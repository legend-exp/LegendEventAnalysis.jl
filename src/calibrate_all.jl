# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).


"""
    calibrate_all(data::LegendData, sel::ValiditySelection, datastore::AbstractDict)

Calibrate all channels in the given datastore, using the metadata
processing configuration for `data` and `sel`.
"""
function calibrate_all(data::LegendData, sel::ValiditySelection, datastore::AbstractDict)
    ds = datastore

    chinfo = channelinfo(data, sel)
    geds_channels = ChannelId.(filterby(@pf $system == :geds && $processable && $usability != :off)(chinfo).channel)
    spms_channels = ChannelId.(filterby(@pf $system == :spms && $processable)(chinfo).channel)
    puls_channels = ChannelId.(filterby(@pf $detector == DetectorId(:PULS01ANA))(chinfo).channel)

    ged_caldata = Dict([
        let detector = channelinfo(data, sel, channel).detector,
            chdata = ds[channel][:]
            channel => _calibrate_ged_data(data, sel, detector, chdata)
        end
        for channel in geds_channels
    ])

    ged_events_pre = build_global_events(ged_caldata, geds_channels)

    # ToDo: Make isvalid_t0 configurable:
    function min_t0(t0::AbstractVector{<:Number})
        isvalid_t0(t0) = 47000u"ns" < t0 < 55000u"ns"
        invalid_to_inf(t0::Number) = !isvalid_t0(t0) ? typeof(t0)(Inf) : t0
        mt = minimum(invalid_to_inf, t0)
        isinf(mt) ? typeof(mt)(NaN) : mt
    end

    emax_ch = broadcast(ged_events_pre.e_trap_cal) do e_cal
        findmax(x -> isnan(x) ? zero(x) : x, e_cal)[2]
    end

    ged_additional_cols = (
        t0_start = [min_t0(t0) for t0 in ged_events_pre.t0],
        multiplicity = [count(e -> (isnan(e) ? zero(e) : e) >= 25u"keV", E) for E in ged_events_pre.e_cusp_cal],
        emax_ch = emax_ch,
        emax_isgood = getindex.(ged_events_pre.isgood, emax_ch),
        emax_trap_cal = getindex.(ged_events_pre.e_trap_cal, emax_ch),
        emax_cusp_cal = getindex.(ged_events_pre.e_trap_cal, emax_ch),
        emax_zac_cal = getindex.(ged_events_pre.e_trap_cal, emax_ch),
        emax_trap_ctc_cal = getindex.(ged_events_pre.e_trap_ctc_cal, emax_ch),
        emax_cusp_ctc_cal = getindex.(ged_events_pre.e_cusp_ctc_cal, emax_ch),
        emax_zac_ctc_cal = getindex.(ged_events_pre.e_zac_ctc_cal, emax_ch),
    )
    ged_events = StructVector(merge(columns(ged_events_pre), ged_additional_cols))
    
    lar_caldata = Dict([
        let detector = channelinfo(data, sel, channel).detector,
            chdata = ds[channel][:]
            channel => _calibrate_spm_data(data, sel, detector, chdata)
        end
        for channel in spms_channels
    ])
    lar_events_novov = build_global_events(lar_caldata, spms_channels)
    lar_events = StructArray(map(_fix_vov, columns(lar_events_novov)))
    
    system_events = (
        geds = ged_events,
        spms = lar_events,
        puls = build_global_events(ds, puls_channels),
    )

    global_events = build_cross_system_events(system_events)

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
