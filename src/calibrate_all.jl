# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).


"""
    calibrate_all(data::LegendData, sel::ValiditySelection, datastore::AbstractDict)

Calibrate all channels in the given datastore, using the metadata
processing configuration for `data` and `sel`.
"""
function calibrate_all(data::LegendData, sel::AnyValiditySelection, datastore::AbstractDict)
    ds = datastore

    chinfo = channelinfo(data, sel)
    geds_channels = ChannelId.(filterby(@pf $system == :geds && $processable && $usability != :off)(chinfo).channel)
    spms_channels = ChannelId.(filterby(@pf $system == :spms && $processable)(chinfo).channel)
    puls_channels = ChannelId.(filterby(@pf $detector == DetectorId(:PULS01ANA))(chinfo).channel)


    # HPGe:

    ged_caldata = Dict([
        let detector = channelinfo(data, sel, channel).detector,
            chdata = ds[channel][:]
            channel => calibrate_ged_channel_data(data, sel, detector, chdata)
        end
        for channel in geds_channels
    ])

    ged_events_pre = build_global_events(ged_caldata, geds_channels)

    min_t0(t0::AbstractVector{<:Number}) = isempty(t0) ? eltype(t0)(NaN) : minimum(t0)

    max_e_ch = broadcast(ged_events_pre.e_cusp_ctc_cal) do e_cal
        findmax(x -> isnan(x) ? zero(x) : x, e_cal)[2]
    end

    trig_e_ch = findall.(ged_events_pre.is_physical_trig)

    trig_e_trap_cal = _fix_vov(getindex.(ged_events_pre.e_trap_cal, trig_e_ch))
    trig_e_cusp_cal = _fix_vov(getindex.(ged_events_pre.e_cusp_cal, trig_e_ch))
    trig_e_zac_cal = _fix_vov(getindex.(ged_events_pre.e_zac_cal, trig_e_ch))
    trig_e_trap_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_trap_ctc_cal, trig_e_ch))
    trig_e_cusp_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_cusp_ctc_cal, trig_e_ch))
    trig_e_zac_ctc_cal = _fix_vov(getindex.(ged_events_pre.e_zac_ctc_cal, trig_e_ch))
    trig_e_short_cal = _fix_vov(getindex.(ged_events_pre.e_short_cal, trig_e_ch))
    trig_t0 = _fix_vov(getindex.(ged_events_pre.t0, trig_e_ch))
    n_trig = length.(trig_e_ch)
    n_expected_baseline = Ref(length(first(ged_events_pre).is_baseline)) .- count.(ged_events_pre.is_baseline)
    

    ged_additional_cols = (
        t0_start = min_t0.(trig_t0),
        trig_t0 = trig_t0,
        multiplicity = n_trig,
        max_e_ch = max_e_ch,
        max_e_trap_cal = maximum.(trig_e_trap_cal, init=0u"keV"),
        max_e_cusp_cal = maximum.(trig_e_cusp_cal, init=0u"keV"),
        max_e_zac_cal = maximum.(trig_e_zac_cal, init=0u"keV"),
        max_e_trap_ctc_cal = maximum.(trig_e_trap_ctc_cal, init=0u"keV"),
        max_e_cusp_ctc_cal = maximum.(trig_e_cusp_ctc_cal, init=0u"keV"),
        max_e_zac_ctc_cal = maximum.(trig_e_zac_ctc_cal, init=0u"keV"),
        max_e_short_cal = maximum.(trig_e_short_cal, init=0u"keV"),
        trig_e_ch = trig_e_ch,
        trig_e_trap_cal = trig_e_trap_cal,
        trig_e_cusp_cal = trig_e_cusp_cal,
        trig_e_zac_cal = trig_e_zac_cal,
        trig_e_trap_ctc_cal = trig_e_trap_ctc_cal,
        trig_e_cusp_ctc_cal = trig_e_cusp_ctc_cal,
        trig_e_zac_ctc_cal = trig_e_zac_ctc_cal,
        trig_e_short_cal = trig_e_short_cal,
        is_valid_qc = count.(ged_events_pre.is_baseline) .== n_expected_baseline .&& count.(ged_events_pre.is_physical) .== n_trig,
    )
    ged_events = StructVector(merge(columns(ged_events_pre), ged_additional_cols))


    # SiPM:

    spm_caldata = Dict([
        let detector = channelinfo(data, sel, channel).detector,
            chdata = ds[channel][:]
            channel => calibrate_spm_channel_data(data, sel, detector, chdata)
        end
        for channel in spms_channels
    ])
    spm_events_novov = build_global_events(spm_caldata, spms_channels)
    spm_events = StructArray(map(_fix_vov, columns(spm_events_novov)))


    # Pulser:

    pls_caldata = Dict([
        let detector = channelinfo(data, sel, channel).detector,
            chdata = ds[channel][:]
            channel => calibrate_pls_channel_data(data, sel, detector, chdata)
        end
        for channel in puls_channels
    ])
    pls_events = build_global_events(pls_caldata, puls_channels)
    #pls_events = map(Broadcast.BroadcastFunction(only), columns(pls_events_pre))
    

    # Cross-system:

    system_events = (
        geds = ged_events,
        spms = spm_events,
        puls = pls_events,
    )

    global_events_pre = build_cross_system_events(system_events)
    single_pls_col = StructVector(map(Broadcast.BroadcastFunction(only), columns(global_events_pre.puls)))
    global_events = StructVector(merge(columns(global_events_pre), (puls = single_pls_col,)))

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
