# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_spm_channel_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, channel_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for the given SiPM
`detector` to the single-channel `channel_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_spm_channel_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, channel_data::AbstractVector; 
    keep_chdata::Bool=false)
    chdata = channel_data[:]

    spmcal_pf = get_spm_cal_propfunc(data, sel, detector)

    # get additional cols to be parsed into the event tier
    chdata_output_pf = if keep_chdata
        PropSelFunction{propertynames(chdata)}()
    else
        get_spms_evt_chdata_propfunc(data, sel)
    end

    cal_output_novv = spmcal_pf.(chdata)
    cal_output = StructArray(map(VectorOfArrays, columns(cal_output_novv)))

    chdata_output = chdata_output_pf.(chdata)

    return StructVector(merge(columns(cal_output), columns(chdata_output)))
end
export calibrate_spm_channel_data


function _single_fiber_esum(
    t_win::AbstractInterval, spm_t::AbstractVector{<:Number},
    spmdc::AbstractVector{<:Number}, spm_pe::AbstractVector{<:Number}
)
    s::eltype(spm_pe) = zero(eltype(spm_pe))
    for i in eachindex(spm_t)
        if spm_t[i] in t_win && !spmdc[i]
            s += spm_pe[i]
        end
    end
    return s
end


# ToDo: Make cut criteria configurable:
function _lar_cut(
    t_win::AbstractInterval, spm_t::AbstractVector{<:AbstractVector{<:Number}},
    spmdc::AbstractVector{<:AbstractVector{<:Number}}, spm_pe::AbstractVector{<:AbstractVector{<:Number}}
)
    n_over_thresh::Int = 0
    pe_sum::eltype(eltype(spm_pe)) = zero(eltype(eltype(spm_pe)))
    for i in eachindex(spm_t)
        s_i = _single_fiber_esum(t_win, spm_t[i], spmdc[i], spm_pe[i])
        if s_i > 0.5
            n_over_thresh += 1
        end
        pe_sum += s_i
    end

    lar_cut = n_over_thresh >= 4 || pe_sum >= 4

    return (lar_cut = lar_cut, spms_win_pe_sum, spms_win_multiplicity = n_over_thresh)
end


# ToDo: Make time window configurable:
function _build_lar_cut(global_events::AbstractVector{<:NamedTuple})
    geds_t0 = global_events.geds.t0_start
    spm_t = global_events.spms.trig_pos
    spmdc = global_events.spms.trig_is_dc
    spm_pe = global_events.spms.trig_pe
    
    t_wins = @. ClosedInterval(geds_t0 - 1u"μs", geds_t0 + 5u"μs")
    return StructArray(_lar_cut.(t_wins, spm_t, spmdc, spm_pe))
end
