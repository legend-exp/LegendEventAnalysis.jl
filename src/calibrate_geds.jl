# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_ged_detector_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, detector_data::AbstractVector; 
        psd_cal_pars_type::Symbol=:ppars, psd_cal_pars_cat::Symbol=:aoe, psd_cut_pars_type::Symbol=:ppars, psd_cut_pars_cat::Symbol=:aoe,
        keep_detdata::Bool=false)

Apply the calibration specified by `data` and `sel` for the given HPGe
`detector` to the single-detector `detector_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_ged_detector_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorIdLike, detector_data::AbstractVector; 
        e_cal_pars_type::Symbol=:rpars, e_cal_pars_cat::Symbol=:ecal,
        psd_pars_type::Symbol=:ppars,
        aoe_cal_pars_type::Symbol=:ppars, aoe_cal_pars_cat::Symbol=:aoe, aoe_cut_pars_type::Symbol=:ppars, aoe_cut_pars_cat::Symbol=:aoe,
        lq_cal_pars_type::Symbol=:ppars, lq_cal_pars_cat::Symbol=:lq, lq_cut_pars_type::Symbol=:ppars, lq_cut_pars_cat::Symbol=:lq,
        keep_detdata::Bool=false)
    
    detector = DetectorId(detector)

    # get all detdata
    detdata = detector_data[:]

    # get energy and psd calibration functions for the detector
    cal_pf = get_ged_cal_propfunc(data, sel, detector; pars_type=e_cal_pars_type, pars_cat=e_cal_pars_cat)
    psd_pf = get_ged_psd_propfunc(data, sel, detector; aoe_pars_type=aoe_cal_pars_type, aoe_pars_cat=aoe_cal_pars_cat, lq_pars_type=lq_cal_pars_type, lq_pars_cat=lq_cal_pars_cat)

    # get qc labels
    cut_pf = get_ged_qc_cuts_propfunc(data, sel, detector)
    
    # get qc cut functions
    cut_is_physical_pf = get_ged_qc_is_physical_propfunc(data, sel, detector)
    cut_is_baseline_pf = get_ged_qc_is_baseline_propfunc(data, sel, detector)
    cut_is_trig_pf = get_ged_qc_is_trig_propfunc(data, sel, detector)
    
    # get additional cols to be parsed into the event tier
    detdata_output_pf = if keep_detdata
        PropSelFunction{propertynames(detdata)}()
    else
        get_ged_evt_detdata_propfunc(data, sel)
    end
    
    # get postcal psd flags
    aoecut_pf = get_ged_aoe_cut_propfunc(data, sel, detector; pars_type=aoe_cut_pars_type, pars_cat=aoe_cut_pars_cat)
    lqcut_pf = get_ged_lq_cut_propfunc(data, sel, detector; pars_type=lq_cut_pars_type, pars_cat=lq_cut_pars_cat)
    psdcut_pf = get_ged_psd_classifier_propfunc(data, sel, detector; pars_type=psd_pars_type)

    # apply calibrations 
    cal_output = cal_pf.(detdata)
    psd_output = psd_pf.(StructArray(merge(columns(cal_output), columns(detdata))))
    cal_detdata = StructArray(merge(columns(cal_output), columns(psd_output), columns(detdata)))

    # get cut labels
    cut_output = cut_pf.(cal_detdata)

    # get postcal data
    postcal_data = StructArray(merge(columns(aoecut_pf.(cal_detdata)), columns(lqcut_pf.(cal_detdata))))
    psd_classifier = psdcut_pf.(postcal_data)

    additional_psd_cols = (
        psd_classifier = psd_classifier,
    )

    # get additional columns 
    detdata_output = detdata_output_pf.(detdata)

    # get qc cut flags
    is_physical = cut_is_physical_pf.(cut_output)
    is_baseline = cut_is_baseline_pf.(cut_output)
    is_physical_trig = cut_is_trig_pf.(cal_detdata) .&& is_physical

    
    additional_qc_cols = (
        is_physical = is_physical,
        is_baseline = is_baseline,
        is_physical_trig = is_physical_trig,
    )

    return StructVector(merge(columns(detdata_output), columns(cal_output), columns(psd_output), columns(postcal_data), additional_psd_cols, columns(cut_output), additional_qc_cols))
end
export calibrate_ged_detector_data

calibrate_ged_detector_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorIdLike, tier::DataTierLike; kwargs...) =
    calibrate_ged_detector_data(data, sel, detector, read_ldata(data, tier, sel, detector), kwargs...)