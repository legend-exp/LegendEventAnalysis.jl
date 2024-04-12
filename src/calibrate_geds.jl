# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_ged_channel_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, channel_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for the given HPGe
`detector` to the single-channel `channel_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_ged_channel_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, channel_data::AbstractVector)
    chdata = channel_data[:]

    # get energy and psd calibration functions for the detector
    cal_pf = get_ged_cal_propfunc(data, sel, detector)

    # get qc labels
    cut_pf = get_ged_qc_cuts_propfunc(data, sel, detector)
    
    # get qc cut functions
    cut_is_physical_pf = get_ged_qc_is_physical_propfunc(data, sel, detector)
    cut_is_baseline_pf = get_ged_qc_is_baseline_propfunc(data, sel, detector)
    cut_is_trig_pf = get_ged_qc_is_trig_propfunc(data, sel, detector)
    
    # get additional cols to be parsed into the event tier
    chdata_output_pf = get_ged_evt_chdata_propfunc(data, sel, detector)
    
    # get postcal psd flags
    postcal_pf = let aoe_window = dataprod_pars_aoe_window(data, sel, detector)
        @pf (
            aeo_low_cut = $aoe_classifier > leftendpoint(aoe_window),
            aoe_ds_cut = $aoe_classifier in aoe_window,
        )
    end
    
    # apply calibrations 
    cal_output = cal_pf.(chdata)
    cal_chdata = StructArray(merge(columns(cal_output), columns(chdata)))

    # get cut labels
    cut_output = cut_pf.(cal_chdata)

    # get additional columns
    chdata_output = chdata_output_pf.(chdata)

    # get postcal data
    postcal_data = postcal_pf.(cal_chdata)

    # get cut flags
    is_physical = cut_is_physical_pf.(cut_output)
    is_baseline = cut_is_baseline_pf.(cut_output)
    is_physical_trig = cut_is_trig_pf.(cal_chdata) .&& is_physical
    
    additional_cols = (
        is_physical = is_physical,
        is_baseline = is_baseline,
        is_physical_trig = is_physical_trig,
    )

    return StructVector(merge(columns(chdata_output), columns(cal_output), columns(postcal_data), columns(cut_output), additional_cols))
end
export calibrate_ged_channel_data
