# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_aux_detector_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, detector_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for aux detector referred to
by the `detector` ID to the single-detector `detector_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_aux_detector_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, detector_data::AbstractVector)
    detdata = detector_data[:]

    # get calibration function for the detector
    auxcal_pf = get_aux_cal_propfunc(data, sel, detector)

    # get additional cols to be parsed into the event tier
    detdata_output_pf = get_aux_evt_detdata_propfunc(data, sel, detector)

    # get additional columns
    detdata_output = detdata_output_pf.(detdata)

    # apply calibrations
    cal_output = auxcal_pf.(detdata)

    return StructVector(merge(columns(detdata_output), columns(cal_output)))
end
export calibrate_aux_detector_data
