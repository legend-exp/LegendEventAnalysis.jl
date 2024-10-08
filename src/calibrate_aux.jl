# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_aux_channel_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, channel_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for aux channel referred to
by the `detector` ID to the single-channel `channel_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_aux_channel_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, channel_data::AbstractVector)
    chdata = channel_data[:]

    # get calibration function for the detector
    auxcal_pf = get_aux_cal_propfunc(data, sel, detector)

    # get additional cols to be parsed into the event tier
    chdata_output_pf = get_aux_evt_chdata_propfunc(data, sel, detector)

    # get additional columns
    chdata_output = chdata_output_pf.(chdata)

    # apply calibrations
    cal_output = auxcal_pf.(chdata)

    return StructVector(merge(columns(chdata_output), columns(cal_output)))
end
export calibrate_aux_channel_data
