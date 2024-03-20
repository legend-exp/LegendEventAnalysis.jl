# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_ged_channel_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, channel_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for the given HPGe
`detector` to the single-channel `channel_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_ged_channel_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, channel_data::AbstractVector)
    chdata = channel_data[:]

    cal_pf = get_ged_cal_propfunc(data, sel, detector)

    # ToDo: Switch to non-hardcoded cuts, see above:
    cut_pf = LegendDataManagement.get_ged_qc_cuts_propfunc(data, sel)
    cut_is_physical_pf = LegendDataManagement.get_ged_qc_is_physical_propfunc(data, sel)
    cut_is_baseline_pf = LegendDataManagement.get_ged_qc_is_baseline_propfunc(data, sel)
    cut_is_trig_pf = LegendDataManagement.get_ged_qc_is_trig_propfunc(data, sel)
    
    # ToDo: Make channel output configurable:
    chdata_output_pf = @pf (
        timestamp = $timestamp,
        t0 = $t0,
        drift_time = $drift_time,
        e_10410 = $e_10410,
        blsigma = $blsigma,
        blmean = $blmean,
    )
    
    postcal_pf = let aoe_window = LegendDataManagement.dataprod_pars_aoe_window(data, sel, detector)
        @pf (
            aeo_low_cut = $aoe_classifier > leftendpoint(aoe_window),
            aoe_ds_cut = $aoe_classifier in aoe_window,
        )
    end
    
    cal_output = cal_pf.(chdata)
    cal_chdata = StructArray(merge(columns(cal_output), columns(chdata)))
    cut_output = cut_pf.(cal_chdata)
    chdata_output = chdata_output_pf.(chdata)
    postcal_data = postcal_pf.(cal_chdata)
    is_physical = cut_is_physical_pf.(cut_output)

    additional_cols = (
        is_physical = is_physical,
        is_baseline = cut_is_baseline_pf.(cut_output),
        is_physical_trig = cut_is_trig_pf.(cal_chdata) .&& is_physical,
    )

    return StructVector(merge(columns(chdata_output), columns(cal_output), columns(postcal_data), columns(cut_output), additional_cols))
end
export calibrate_ged_channel_data
