# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).


"""
    calibrate_all(data::LegendData, sel::ValiditySelection, datastore::AbstractDict)

Calibrate all channels in the given datastore, using the metadata
processing configuration for `data` and `sel`.
"""
function calibrate_all(data::LegendData, sel::AnyValiditySelection, datastore::AbstractDict, tier::DataTierLike=:jldsp)
    ds = datastore

    @debug "Calibrating all channels in for `ValiditySelection` $(sel) in `DataTier` $(tier)"
    chinfo = channelinfo(data, sel)
    geds_channels::Vector{ChannelId} = filterby(get_ged_evt_chsel_propfunc(data, sel))(chinfo).channel
    @debug "Loaded $(length(geds_channels)) HPGe channels"
    spms_channels::Vector{ChannelId} = filterby(get_spms_evt_chsel_propfunc(data, sel))(chinfo).channel
    @debug "Loaded $(length(spms_channels)) SiPM channels"
    aux_channels::Vector{ChannelId} = filterby(get_aux_evt_chsel_propfunc(data, sel))(chinfo).channel
    @debug "Loaded auxiliary channels: $(join(string.(filterby(get_aux_evt_chsel_propfunc(data, sel))(chinfo).detector), ", "))"

    # Main.@infiltrate
    # HPGe:
    @debug "Calibrating HPGe channels"
    ged_kwargs = get_ged_evt_kwargs(data, sel)
    ged_caldata_v = Vector{StructVector}(undef, length(geds_channels))
    Threads.@threads for i in eachindex(geds_channels)
        let detector = channelinfo(data, sel, geds_channels[i]).detector, chdata = ds[geds_channels[i], tier][:]
            ged_caldata_v[i] = calibrate_ged_channel_data(data, sel, detector, chdata; ged_kwargs...)
        end
    end
    ged_caldata = Dict(geds_channels .=> ged_caldata_v)

    @debug "Building global events for HPGe channels"
    ged_events_pre = build_global_events(ged_caldata, geds_channels)

    min_t0(t0::AbstractVector{<:Number}) = isempty(t0) ? eltype(t0)(NaN) : minimum(t0)

    max_e_ch = broadcast(ged_events_pre.e_cusp_ctc_cal) do e_cal
        findmax(x -> isnan(x) ? zero(x) : x, e_cal)[2]
    end

    trig_e_ch = findall.(ged_events_pre.is_physical_trig)

    trig_e_trap_cal = _fix_vov(getindex.(ged_events_pre.e_trap_cal, trig_e_ch))
    trig_e_cusp_cal = _fix_vov(getindex.(ged_events_pre.e_cusp_cal, trig_e_ch))
    trig_e_trap_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_trap_ctc_cal, trig_e_ch))
    trig_e_cusp_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_cusp_ctc_cal, trig_e_ch))
    trig_e_short_cal = _fix_vov(getindex.(ged_events_pre.e_313_cal, trig_e_ch))
    trig_t0 = _fix_vov(getindex.(ged_events_pre.t0, trig_e_ch))
    n_trig = length.(trig_e_ch)
    n_expected_baseline = length.(ged_events_pre.is_baseline) .- length.(trig_e_ch)
    
    maximum_with_init(A) = maximum(A, init=zero(eltype((A))))

    ged_additional_cols = (
        t0_start = min_t0.(trig_t0),
        trig_t0 = trig_t0,
        multiplicity = n_trig,
        max_e_ch = max_e_ch,
        max_e_trap_cal = maximum_with_init.(trig_e_trap_cal),
        max_e_cusp_cal = maximum_with_init.(trig_e_cusp_cal),
        max_e_trap_ctc_cal = maximum_with_init.(trig_e_trap_ctc_cal),
        max_e_cusp_ctc_cal = maximum_with_init.(trig_e_cusp_ctc_cal),
        max_e_short_cal = maximum_with_init.(trig_e_short_cal),
        trig_e_ch = trig_e_ch,
        trig_e_trap_cal = trig_e_trap_cal,
        trig_e_cusp_cal = trig_e_cusp_cal,
        trig_e_trap_ctc_cal = trig_e_trap_ctc_cal,
        trig_e_cusp_ctc_cal = trig_e_cusp_ctc_cal,
        trig_e_short_cal = trig_e_short_cal,
        is_valid_qc = count.(ged_events_pre.is_baseline) .== n_expected_baseline,
        is_discharge_recovery = any.(ged_events_pre.is_discharge_recovery_ml),
        is_saturated = any.(ged_events_pre.is_saturated),
        is_discharge = any.(ged_events_pre.is_discharge),
    )
    ged_events = StructVector(merge(columns(ged_events_pre), ged_additional_cols))


    # SiPM:
    @debug "Calibrating SiPM channels"
    spm_caldata_v = Vector{StructVector}(undef, length(spms_channels))
    Threads.@threads for i in eachindex(spms_channels)
        let detector = channelinfo(data, sel, spms_channels[i]).detector, chdata = ds[spms_channels[i], tier][:]
            spm_caldata_v[i] = calibrate_spm_channel_data(data, sel, detector, chdata)
        end
    end
    spm_caldata = Dict(spms_channels .=> spm_caldata_v)
    @debug "Building global events for SiPM channels"
    spm_events_novov = build_global_events(spm_caldata, spms_channels)
    spm_events = StructArray(map(_fix_vov, columns(spm_events_novov)))

    
    @debug "Calibrating auxiliary channels"
    # aux & Forced Trigger
    aux_caldata = 
        [Dict(
            let detector = channelinfo(data, sel, channel).detector,
                chdata = ds[channel, tier][:]
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
        ged_spm = _build_lar_cut(global_events),
    )

    result = StructArray(merge(columns(global_events), cross_systems_cols))
    
    return result
end
export calibrate_all


_fix_vov(x) = x
_fix_vov(x::AbstractVector{<:AbstractVector}) = VectorOfVectors(x)
_fix_vov(x::VectorOfVectors{<:AbstractVector}) = VectorOfVectors(VectorOfVectors(flatview(x)), x.elem_ptr)
