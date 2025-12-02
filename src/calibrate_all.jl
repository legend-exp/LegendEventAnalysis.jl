# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).


"""
    _get_channel_data(ds::AbstractDataStore, data::LegendData, sel::AnyValiditySelection, channel::ChannelId, tier::DataTierLike)

Helper function for flexible key resolution when accessing channel data from datastore.
Prefers detector name key (e.g., "V00048A"), falls back to channel ID key (e.g., "ch1084000").
This ensures compatibility with both old (channel ID) and new (detector name) file formats.
"""
function _get_channel_data(ds::AbstractDataStore, data::LegendData, sel::AnyValiditySelection, channel::ChannelId, tier::DataTierLike)
    detector = channelinfo(data, sel, channel).detector
    detector_key = string(detector)
    channel_key = string(channel)
    
    # Try detector name first (new format)
    if haskey(ds, detector_key)
        return detector, ds[detector, tier][:]
    # Fallback to channel ID (old format)
    elseif haskey(ds, channel_key)
        return detector, ds[channel, tier][:]
    else
        throw(KeyError("Neither detector key '$detector_key' nor channel key '$channel_key' found in datastore"))
    end
end

"""
    _has_channel_data(ds::AbstractDataStore, data::LegendData, sel::AnyValiditySelection, channel::ChannelId)

Check if channel data exists in datastore under either detector name or channel ID key.
"""
function _has_channel_data(ds::AbstractDataStore, data::LegendData, sel::AnyValiditySelection, channel::ChannelId)
    detector = channelinfo(data, sel, channel).detector
    detector_key = string(detector)
    channel_key = string(channel)
    return haskey(ds, detector_key) || haskey(ds, channel_key)
end


"""
    calibrate_all(data::LegendData, sel::ValiditySelection, datastore::AbstractDict)

Calibrate all channels in the given datastore, using the metadata
processing configuration for `data` and `sel`.
"""
function calibrate_all(data::LegendData, sel::AnyValiditySelection, datastore::AbstractDataStore, tier::DataTierLike=:jldsp)
    ds = datastore

    @debug "Calibrating all channels in for `ValiditySelection` $(sel) in `DataTier` $(tier)"
    chinfo = channelinfo(data, sel)
    geds_channels::Vector{ChannelId} = filterby(get_ged_evt_chsel_propfunc(data, sel))(chinfo).channel
    @debug "Loaded $(length(geds_channels)) HPGe channels"
    hitgeds_channels::Vector{ChannelId} = filterby(get_ged_evt_hitchsel_propfunc(data, sel))(chinfo).channel
    @debug "Loaded $(length(hitgeds_channels)) HPGe hit channels"
    spms_channels::Vector{ChannelId} = filterby(get_spms_evt_chsel_propfunc(data, sel))(chinfo).channel
    @debug "Loaded $(length(spms_channels)) SiPM channels"
    pmts_channels::Vector{ChannelId} = filterby(get_pmts_evt_chsel_propfunc(data, sel))(chinfo).channel
    @debug "Loaded $(length(pmts_channels)) PMT channels"
    aux_channels::Vector{ChannelId} = filterby(get_aux_evt_chsel_propfunc(data, sel))(chinfo).channel
    @debug "Loaded auxiliary channels: $(join(string.(filterby(get_aux_evt_chsel_propfunc(data, sel))(chinfo).detector), ", "))"

    # HPGe:
    @debug "Calibrating HPGe channels"
    ged_kwargs = get_ged_evt_kwargs(data, sel)
    ged_caldata_v = Vector{StructVector}(undef, length(geds_channels))
    p = Progress(length(geds_channels); desc="Calibrating HPGe channels...")
    Threads.@threads for i in eachindex(geds_channels)
        let (detector, chdata) = _get_channel_data(ds, data, sel, geds_channels[i], tier)
            ged_caldata_v[i] = calibrate_ged_channel_data(data, sel, detector, chdata; ged_kwargs...)
            next!(p; showvalues = [("Calibrated detector", detector)])
        end
    end
    ged_caldata = Dict(geds_channels .=> ged_caldata_v)

    @debug "Building global events for HPGe channels"
    ged_events_pre = build_global_events(ged_caldata, geds_channels)

    # Get optional t0 valid window from PMT muon_cut config (filters out bad t0 values from AC-coupled detectors)
    dataprod_muoncut = get_pmts_evt_muon_cut_props(data, sel)
    ged_t0_valid_window = get(dataprod_muoncut, :ged_t0_valid_window, nothing)
    # Note: ged_t0_valid_window is a tuple of Quantity values like (40.0µs, 60.0µs)
    # t0 values from DSP also have units (µs), so comparison works directly
    
    # Filter out NaN values and optionally values outside valid window when computing minimum t0
    function min_t0_filtered(t0::AbstractVector, ged_t0_window)
        # Filter NaN values (works with both unitful and unitless)
        valid_t0 = filter(x -> !isnan(ustrip(x)), t0)
        if !isnothing(ged_t0_window)
            t0_min, t0_max = first(ged_t0_window), last(ged_t0_window)
            # Strip units for comparison to handle both unitful and unitless t0 values
            t0_min_val = ustrip(u"µs", t0_min)
            t0_max_val = ustrip(u"µs", t0_max)
            valid_t0 = filter(t -> ustrip(u"µs", t) >= t0_min_val && ustrip(u"µs", t) <= t0_max_val, valid_t0)
        end
        return isempty(valid_t0) ? eltype(t0)(NaN) : minimum(valid_t0)
    end

    max_e_ch = broadcast(ged_events_pre.e_cusp_ctc_cal) do e_cal
        findmax(x -> isnan(x) ? zero(x) : x, e_cal)[2]
    end

    trig_e_ch = findall.(ged_events_pre.is_physical_trig)

    trig_e_trap_max_cal = _fix_vov(getindex.(ged_events_pre.e_trap_max_cal, trig_e_ch))
    trig_e_cusp_max_cal = _fix_vov(getindex.(ged_events_pre.e_cusp_max_cal, trig_e_ch))
    trig_e_trap_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_trap_ctc_cal, trig_e_ch))
    trig_e_cusp_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_cusp_ctc_cal, trig_e_ch))
    trig_e_535_cal      = _fix_vov(getindex.(ged_events_pre.e_535_cal, trig_e_ch))
    trig_t0 = _fix_vov(getindex.(ged_events_pre.t0, trig_e_ch))

    # Check if any events have only NaN t0 values (detectors without calibration)
    # Note: Use ustrip for Unitful compatibility
    n_events_all_nan_t0 = count(t0 -> all(x -> isnan(ustrip(x)), t0), trig_t0)
    if n_events_all_nan_t0 > 0
        @warn "$(n_events_all_nan_t0) events have all NaN t0 values (all trigger detectors missing calibration). " *
              "t0_start will be NaN for these events, affecting muon veto coincidence matching."
    end

    n_trig = length.(trig_e_ch)
    n_expected_baseline = length.(ged_events_pre.is_baseline) .- length.(trig_e_ch)
    
    maximum_with_init(A) = maximum(A, init=zero(eltype((A))))

    is_valid_trig(trig_chs::AbstractVector{<:Int}, hit_channels::AbstractVector{<:Int}) = all(x -> x in hit_channels, trig_chs)
    
    hit_properties = get_ged_evt_is_valid_hit_properties(data, sel)
    is_valid_hit = trues(length(ged_events_pre.channel))
    for prop in hit_properties
        is_valid_hit .&= all.(map.(isfinite, (getindex.(getproperty(ged_events_pre, prop), trig_e_ch))))
    end

    ged_additional_cols = (
        t0_start = min_t0_filtered.(trig_t0, Ref(ged_t0_valid_window)),
        trig_t0 = trig_t0,
        multiplicity = n_trig,
        max_e_ch_idxs = max_e_ch,
        max_e_ch = only.(getindex.(ged_events_pre.channel, max_e_ch)),
        max_e_trap_cal = maximum_with_init.(trig_e_trap_max_cal),
        max_e_cusp_cal = maximum_with_init.(trig_e_cusp_max_cal),
        max_e_trap_ctc_cal = maximum_with_init.(trig_e_trap_ctc_cal),
        max_e_cusp_ctc_cal = maximum_with_init.(trig_e_cusp_ctc_cal),
        max_e_short_cal = maximum_with_init.(trig_e_535_cal),
        trig_e_ch_idxs = trig_e_ch,
        trig_e_ch = getindex.(ged_events_pre.channel, trig_e_ch),
        trig_e_trap_max_cal = trig_e_trap_max_cal,
        trig_e_cusp_max_cal = trig_e_cusp_max_cal,
        trig_e_trap_ctc_cal = trig_e_trap_ctc_cal,
        trig_e_cusp_ctc_cal = trig_e_cusp_ctc_cal,
        trig_e_535_cal = trig_e_535_cal,
        is_valid_qc = count.(ged_events_pre.is_baseline) .== n_expected_baseline,
        is_valid_trig = is_valid_trig.(getindex.(ged_events_pre.channel, trig_e_ch), Ref(Int.(hitgeds_channels))),
        is_valid_hit = is_valid_hit,
        is_valid_psd = all.(getindex.(ged_events_pre.psd_classifier, trig_e_ch)),
        is_discharge_recovery = any.(ged_events_pre.is_discharge_recovery_ml),
        is_saturated = any.(ged_events_pre.is_saturated),
        is_discharge = any.(ged_events_pre.is_discharge),
    )
    ged_events = StructVector(merge(columns(ged_events_pre), ged_additional_cols))


    # SiPM:
    @debug "Calibrating SiPM channels"
    spm_kwargs = get_spms_evt_kwargs(data, sel)
    spm_caldata_v = Vector{StructVector}(undef, length(spms_channels))
    p = Progress(length(spms_channels); desc="Calibrating SiPM channels...")
    Threads.@threads for i in eachindex(spms_channels)
        let (detector, chdata) = _get_channel_data(ds, data, sel, spms_channels[i], tier)
            spm_caldata_v[i] = calibrate_spm_channel_data(data, sel, detector, chdata; spm_kwargs...)
            next!(p; showvalues = [("Calibrated detector", detector)])
        end
    end
    spm_caldata = Dict(spms_channels .=> spm_caldata_v)
    @debug "Building global events for SiPM channels"
    spm_events_novov = build_global_events(spm_caldata, spms_channels)
    spm_events = StructArray(map(_fix_vov, columns(spm_events_novov)))

    # PMT:
    pmt_caldata = nothing  # Will be set if PMT data exists
    pmt_events = if all(.!_has_channel_data.(Ref(ds), Ref(data), Ref(sel), pmts_channels))
        @warn "No PMT data found, skip PMT calibration"
        Vector{NamedTuple{(:timestamp, ), Tuple{Unitful.Time{<:Real}, }}}()
    else
        @debug "Calibrating PMT channels"
        pmt_kwargs = get_pmts_evt_kwargs(data, sel)
        pmt_caldata_v = Vector{StructVector}(undef, length(pmts_channels))
        p = Progress(length(pmts_channels); desc="Calibrating PMT channels...")
        Threads.@threads for i in eachindex(pmts_channels)
            let (detector, chdata) = _get_channel_data(ds, data, sel, pmts_channels[i], tier)
                pmt_caldata_v[i] = calibrate_pmt_channel_data(data, sel, detector, chdata; pmt_kwargs...)
                next!(p; showvalues = [("Calibrated detector", detector)])
            end
        end
        pmt_caldata = Dict(pmts_channels .=> pmt_caldata_v)
        @debug "Building global events for PMT channels"
        pmt_events_pre_novov = build_global_events(pmt_caldata, pmts_channels)
        pmt_events_pre = StructArray(map(_fix_vov, columns(pmt_events_pre_novov)))

        StructVector(merge(columns(_build_muon_cut(data, sel, pmt_events_pre)), columns(pmt_events_pre)))
    end


    @debug "Calibrating auxiliary channels"
    # aux & Forced Trigger
    aux_caldata = 
        [Dict(
            let (detector, chdata) = _get_channel_data(ds, data, sel, channel, tier)
                channel => calibrate_aux_channel_data(data, sel, detector, chdata)
            end
            ) for channel in aux_channels]
    @debug "Building global events for auxiliary channels"
    aux_events = NamedTuple{Tuple(get_aux_evt_levelname_propfunc.(Ref(data), Ref(sel), reduce(vcat, collect.(keys.(aux_caldata)))))}(build_global_events.(aux_caldata))

    # Cross-system:
    @debug "Building cross-system events"
    system_events = merge((
        geds = ged_events,
        spms = spm_events,
    ), aux_events)

    global_events_pre = build_cross_system_events(system_events)
    aux_cols = NamedTuple{keys(aux_events)}([StructVector(map(Broadcast.BroadcastFunction(only), columns(getproperty(global_events_pre, k)))) for k in keys(aux_events)])
    global_events = StructVector(merge(Base.structdiff(columns(global_events_pre), NamedTuple{keys(aux_events)}), (aux = StructArray(aux_cols),)))

    cross_systems_cols = (
        ged_spm = _build_lar_cut(data, sel, global_events, global_events.geds.t0_start),
        ft_spm = _build_lar_cut(data, sel, global_events, fill(get_spms_evt_lar_cut_props(data, sel).ft_cut_t0, length(global_events.geds.t0_start))),
        ged_pmt = _build_muon_evt_cut(data, sel, global_events, pmt_events)
    )

    result = StructArray(merge(columns(global_events), cross_systems_cols))
    
    # Check calibration quality and identify detectors with missing parameters
    calibration_issues = _check_calibration_quality(data, sel, ged_caldata, spm_caldata, pmt_caldata, 
                                                     geds_channels, spms_channels, pmts_channels)
    
    result_t = Table(NamedTuple{propertynames(result)}([if c isa StructArray Table(c) else c end for c in columns(result)]))
    return result_t, Table(pmt_events), calibration_issues
end
export calibrate_all


_fix_vov(x) = x
_fix_vov(x::AbstractVector{<:AbstractVector}) = VectorOfVectors(x)
_fix_vov(x::VectorOfVectors{<:AbstractVector}) = VectorOfVectors(VectorOfVectors(flatview(x)), x.elem_ptr)
_fix_vov(t::Table) = Table(NamedTuple{propertynames(t)}([if c isa Table Table(StructArray(map(_fix_vov, columns(c)))) else c end for c in columns(t)]))


"""
    _check_calibration_quality(data, sel, ged_caldata, spm_caldata, pmt_caldata, geds_channels, spms_channels, pmts_channels)

Identify detectors with missing calibration parameters (NaN fallback values).
Returns a NamedTuple with lists of affected detector names per subsystem.
"""
function _check_calibration_quality(data::LegendData, sel::AnyValiditySelection,
                                     ged_caldata::Dict, spm_caldata::Dict, pmt_caldata::Union{Dict, Nothing},
                                     geds_channels::Vector{ChannelId}, 
                                     spms_channels::Vector{ChannelId}, 
                                     pmts_channels::Vector{ChannelId})
    
    missing_ecal = String[]
    missing_sipmcal = String[]
    missing_pmtcal = String[]
    
    # === HPGe: Check for NaN energy calibration per detector ===
    for ch in geds_channels
        try
            chinfo = channelinfo(data, sel, ch)
            detector = chinfo.detector
            usability = chinfo.usability
            caldata = ged_caldata[ch]
            # Check if e_cusp_ctc_cal has any valid (non-NaN) values
            if hasproperty(caldata, :e_cusp_ctc_cal)
                e_vals = caldata.e_cusp_ctc_cal
                if !isempty(e_vals) && all(isnan, e_vals)
                    # Include usability status in detector name
                    push!(missing_ecal, "$(detector) ($(usability))")
                end
            end
        catch; end
    end
    
    # === SiPM: Check for NaN energy calibration per detector ===
    for ch in spms_channels
        try
            detector = channelinfo(data, sel, ch).detector
            caldata = spm_caldata[ch]
            # Check if trig_max_cal (main SiPM energy) has NaN values
            has_nan = false
            if hasproperty(caldata, :trig_max_cal)
                e_vals = caldata.trig_max_cal
                # trig_max_cal is VectorOfArrays, need to flatten
                flat_vals = reduce(vcat, e_vals; init=Float64[])
                if !isempty(flat_vals) && all(isnan, flat_vals)
                    has_nan = true
                end
            end
            if has_nan
                push!(missing_sipmcal, string(detector))
            end
        catch; end
    end
    
    # === PMT: Check for NaN energy calibration per detector ===
    if pmt_caldata !== nothing
        for ch in pmts_channels
            try
                detector = channelinfo(data, sel, ch).detector
                caldata = pmt_caldata[ch]
                # Check if e_fc has NaN values  
                if hasproperty(caldata, :e_fc)
                    e_vals = caldata.e_fc
                    if !isempty(e_vals) && all(isnan, e_vals)
                        push!(missing_pmtcal, string(detector))
                    end
                end
            catch; end
        end
    else
        # All PMTs missing if no pmt_caldata
        for ch in pmts_channels
            try
                detector = channelinfo(data, sel, ch).detector
                push!(missing_pmtcal, string(detector))
            catch; end
        end
    end
    
    # Log warnings if any detectors have missing calibrations
    if !isempty(missing_ecal)
        @warn "HPGe detectors with missing ecal (NaN energy): $(join(missing_ecal, ", "))"
    end
    if !isempty(missing_sipmcal)
        @warn "SiPM detectors with missing sipmcal: $(join(missing_sipmcal, ", "))"
    end
    if !isempty(missing_pmtcal)
        @warn "PMT detectors with missing pmtcal: $(join(missing_pmtcal, ", "))"
    end
    
    return (
        missing_ecal = missing_ecal,
        missing_sipmcal = missing_sipmcal,
        missing_pmtcal = missing_pmtcal
    )
end