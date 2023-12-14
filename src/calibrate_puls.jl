# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_pls_channel_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, channel_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for pulser referred to
by the `detector` ID to the single-channel `channel_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_pls_channel_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, channel_data::AbstractVector)
    chdata = channel_data[:]

    pulsercal_pf = get_pulser_cal_propfunc(data, sel, detector)

    # ToDo: Make channel output configurable:
    chdata_output = (
        timestamp = chdata.timestamp,
        e_10410 = chdata.e_10410,
        t50 = chdata.t50,
    )

    cal_output = pulsercal_pf.(chdata)

    return StructVector(merge(chdata_output, columns(cal_output)))
end
export calibrate_pls_channel_data
