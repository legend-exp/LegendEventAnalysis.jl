# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).


"""
    calibrate_all(data::LegendData, sel::ValiditySelection, datastore::AbstractDict)

Calibrate all detectors in the given datastore, using the metadata
processing configuration for `data` and `sel`.
"""
function calibrate_all(data::LegendData, sel::AnyValiditySelection, datastore::AbstractDataStore, tier::DataTierLike=:jldsp)
    ds = datastore

    @debug "Calibrating all detectors for `ValiditySelection` $(sel) in `DataTier` $(tier)"
    chinfo = channelinfo(data, sel)
    geds_detectors::Vector{DetectorId} = filterby(get_ged_evt_detsel_propfunc(data, sel))(chinfo).detector
    @debug "Loaded $(length(geds_detectors)) HPGe detectors"
    hitgeds_detectors::Vector{DetectorId} = filterby(get_ged_evt_hitdetsel_propfunc(data, sel))(chinfo).detector
    @debug "Loaded $(length(hitgeds_detectors)) HPGe hit detectors"
    spms_detectors::Vector{DetectorId} = filterby(get_spms_evt_detsel_propfunc(data, sel))(chinfo).detector
    @debug "Loaded $(length(spms_detectors)) SiPM detectors"
    pmts_detectors::Vector{DetectorId} = filterby(get_pmts_evt_detsel_propfunc(data, sel))(chinfo).detector
    @debug "Loaded $(length(pmts_detectors)) PMT detectors"
    aux_detectors::Vector{DetectorId} = filterby(get_aux_evt_detsel_propfunc(data, sel))(chinfo).detector
    @debug "Loaded auxiliary detectors: $(join(string.(aux_detectors), ", "))"

    # HPGe:
    @debug "Calibrating HPGe detectors"
    ged_kwargs = get_ged_evt_kwargs(data, sel)
    ged_caldata_v = Vector{StructVector}(undef, length(geds_detectors))
    p = Progress(length(geds_detectors); desc="Calibrating HPGe detectors...")
    Threads.@threads for i in eachindex(geds_detectors)
        let detector = geds_detectors[i], detdata = ds[string(detector), tier][:]
            ged_caldata_v[i] = calibrate_ged_detector_data(data, sel, detector, detdata; ged_kwargs...)
            next!(p; showvalues = [("Calibrated detector", detector)])
        end
    end
    ged_caldata = Dict(geds_detectors .=> ged_caldata_v)

    @debug "Building global events for HPGe detectors"
    ged_events_pre = build_global_events(ged_caldata, geds_detectors)

    min_t0(t0::AbstractVector{<:Number}) = isempty(t0) ? eltype(t0)(NaN) : minimum(t0)

    max_e_det = broadcast(ged_events_pre.e_cusp_ctc_cal) do e_cal
        findmax(x -> isnan(x) ? zero(x) : x, e_cal)[2]
    end

    trig_e_det = findall.(ged_events_pre.is_physical_trig)

    trig_e_trap_max_cal = _fix_vov(getindex.(ged_events_pre.e_trap_max_cal, trig_e_det))
    trig_e_cusp_max_cal = _fix_vov(getindex.(ged_events_pre.e_cusp_max_cal, trig_e_det))
    trig_e_trap_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_trap_ctc_cal, trig_e_det))
    trig_e_cusp_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_cusp_ctc_cal, trig_e_det))
    trig_e_535_cal      = _fix_vov(getindex.(ged_events_pre.e_535_cal, trig_e_det))
    trig_t0 = _fix_vov(getindex.(ged_events_pre.t0, trig_e_det))
    n_trig = length.(trig_e_det)
    n_expected_baseline = length.(ged_events_pre.is_baseline) .- length.(trig_e_det)
    
    maximum_with_init(A) = maximum(A, init=zero(eltype((A))))

    is_valid_trig(trig_dets::AbstractVector{<:DetectorId}, hit_detectors::AbstractVector{<:DetectorId}) = all(x -> x in hit_detectors, trig_dets)
    
    hit_properties = get_ged_evt_is_valid_hit_properties(data, sel)
    is_valid_hit = trues(length(ged_events_pre.detector))
    for prop in hit_properties
        is_valid_hit .&= all.(map.(isfinite, (getindex.(getproperty(ged_events_pre, prop), trig_e_det))))
    end

    ged_additional_cols = (
        t0_start = min_t0.(trig_t0),
        trig_t0 = trig_t0,
        multiplicity = n_trig,
        max_e_det_idxs = max_e_det,
        max_e_det = getindex.(ged_events_pre.detector, max_e_det),
        max_e_trap_cal = maximum_with_init.(trig_e_trap_max_cal),
        max_e_cusp_cal = maximum_with_init.(trig_e_cusp_max_cal),
        max_e_trap_ctc_cal = maximum_with_init.(trig_e_trap_ctc_cal),
        max_e_cusp_ctc_cal = maximum_with_init.(trig_e_cusp_ctc_cal),
        max_e_short_cal = maximum_with_init.(trig_e_535_cal),
        trig_e_det_idxs = trig_e_det,
        trig_e_det = getindex.(ged_events_pre.detector, trig_e_det),
        trig_e_trap_max_cal = trig_e_trap_max_cal,
        trig_e_cusp_max_cal = trig_e_cusp_max_cal,
        trig_e_trap_ctc_cal = trig_e_trap_ctc_cal,
        trig_e_cusp_ctc_cal = trig_e_cusp_ctc_cal,
        trig_e_535_cal = trig_e_535_cal,
        is_valid_qc = count.(ged_events_pre.is_baseline) .== n_expected_baseline,
        is_valid_trig = is_valid_trig.(getindex.(ged_events_pre.detector, trig_e_det), Ref(hitgeds_detectors)),
        is_valid_hit = is_valid_hit,
        is_valid_psd = all.(getindex.(ged_events_pre.psd_classifier, trig_e_det)),
        is_discharge_recovery = any.(ged_events_pre.is_discharge_recovery_ml),
        is_saturated = any.(ged_events_pre.is_saturated),
        is_discharge = any.(ged_events_pre.is_discharge),
    )
    ged_events = StructVector(merge(columns(ged_events_pre), ged_additional_cols))


    # SiPM:
    @debug "Calibrating SiPM detectors"
    spm_kwargs = get_spms_evt_kwargs(data, sel)
    spm_caldata_v = Vector{StructVector}(undef, length(spms_detectors))
    p = Progress(length(spms_detectors); desc="Calibrating SiPM detectors...")
    Threads.@threads for i in eachindex(spms_detectors)
        let detector = spms_detectors[i], detdata = ds[string(detector), tier][:]
            spm_caldata_v[i] = calibrate_spm_detector_data(data, sel, detector, detdata; spm_kwargs...)
            next!(p; showvalues = [("Calibrated detector", detector)])
        end
    end
    spm_caldata = Dict(spms_detectors .=> spm_caldata_v)
    @debug "Building global events for SiPM detectors"
    spm_events_novov = build_global_events(spm_caldata, spms_detectors)
    spm_events = StructArray(map(_fix_vov, columns(spm_events_novov)))

    # PMT:
    pmt_events = if all(.!haskey.(Ref(ds), string.(pmts_detectors)))
        @warn "No PMT data found, skip PMT calibration"
        Vector{NamedTuple{(:timestamp, ), Tuple{Unitful.Time{<:Real}, }}}()
    else
        @debug "Calibrating PMT detectors"
        pmt_kwargs = get_pmts_evt_kwargs(data, sel)
        pmt_caldata_v = Vector{StructVector}(undef, length(pmts_detectors))
        p = Progress(length(pmts_detectors); desc="Calibrating PMT detectors...")
        Threads.@threads for i in eachindex(pmts_detectors)
            let detector = pmts_detectors[i], detdata = ds[string(detector), tier][:]
                pmt_caldata_v[i] = calibrate_pmt_detector_data(data, sel, detector, detdata; pmt_kwargs...)
                next!(p; showvalues = [("Calibrated detector", detector)])
            end
        end
        pmt_caldata = Dict(pmts_detectors .=> pmt_caldata_v)
        @debug "Building global events for PMT detectors"
        pmt_events_pre_novov = build_global_events(pmt_caldata, pmts_detectors)
        pmt_events_pre = StructArray(map(_fix_vov, columns(pmt_events_pre_novov)))

        StructVector(merge(columns(_build_muon_cut(data, sel, pmt_events_pre)), columns(pmt_events_pre)))
    end


    @debug "Calibrating auxiliary detectors"
    # aux & Forced Trigger
    aux_caldata = 
        [Dict(
            let detector = aux_detectors[i],
                detdata = ds[string(detector), tier][:]
                detector => calibrate_aux_detector_data(data, sel, detector, detdata)
            end
            ) for i in eachindex(aux_detectors)]
    @debug "Building global events for auxiliary detectors"
    aux_events = NamedTuple{Tuple(get_aux_evt_levelname_propfunc.(Ref(data), Ref(sel), reduce(vcat, collect.(keys.(aux_caldata)))))}(build_global_events.(aux_caldata))

    # Cross-system:
    @debug "Building cross-system events"
    system_events = merge((
        geds = ged_events,
        spms = spm_events,
    ), aux_events)

    global_events_pre = build_cross_system_events(system_events)
    _aux_only(c) = length(c) == 1 ? only(c) : (@warn "Duplicate aux entry, using first"; first(c))
    aux_cols = NamedTuple{keys(aux_events)}([StructVector(map(Broadcast.BroadcastFunction(_aux_only), columns(getproperty(global_events_pre, k)))) for k in keys(aux_events)])
    global_events = StructVector(merge(Base.structdiff(columns(global_events_pre), NamedTuple{keys(aux_events)}), (aux = StructArray(aux_cols),)))

    cross_systems_cols = (
        ged_spm = _build_lar_cut(data, sel, global_events, global_events.geds.t0_start),
        ft_spm = _build_lar_cut(data, sel, global_events, fill(get_spms_evt_lar_cut_props(data, sel).ft_cut_t0, length(global_events.geds.t0_start))),
        ged_pmt = _build_muon_evt_cut(data, sel, global_events, pmt_events)
    )

    result = StructArray(merge(columns(global_events), cross_systems_cols))
    
    result_t = Table(NamedTuple{propertynames(result)}([if c isa StructArray Table(c) else c end for c in columns(result)]))
    return result_t, Table(pmt_events)
end
export calibrate_all


_fix_vov(x) = x
_fix_vov(x::AbstractVector{<:AbstractVector}) = VectorOfVectors(x)
_fix_vov(x::VectorOfVectors{<:AbstractVector}) = VectorOfVectors(VectorOfVectors(flatview(x)), x.elem_ptr)
_fix_vov(t::Table) = Table(NamedTuple{propertynames(t)}([if c isa Table Table(StructArray(map(_fix_vov, columns(c)))) else c end for c in columns(t)]))