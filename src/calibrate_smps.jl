# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_spm_detector_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, detector_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for the given SiPM
`detector` to the single-detector `detector_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_spm_detector_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, detector_data::AbstractVector;
    e_cal_pars_type::Symbol=:ppars, e_cal_pars_cat::Symbol=:sipmcal, dc_cut_pars_type::Symbol=:ppars,
    keep_detdata::Bool=false)
    detdata = detector_data[:]

    spmcal_pf = get_spm_cal_propfunc(data, sel, detector; pars_type=e_cal_pars_type, pars_cat=e_cal_pars_cat)
    spmdc_sel_pf = get_spm_dc_sel_propfunc(data, sel, detector; pars_type=dc_cut_pars_type)
    spmdc_cal_pf = get_spm_dc_cal_propfunc(data, sel, detector; pars_type=dc_cut_pars_type)

    # get additional cols to be parsed into the event tier
    detdata_output_pf = if keep_detdata
        PropSelFunction{propertynames(detdata)}()
    else
        get_spms_evt_detdata_propfunc(data, sel)
    end

    cal_output_novv = spmcal_pf.(detdata)
    cal_output = StructArray(map(VectorOfArrays, columns(cal_output_novv)))

    dc_output_novv = NamedTuple{keys(spmdc_sel_pf)}([spmdc_cal_pf[e_type].(spmdc_sel_pf[e_type].(detdata)) for e_type in keys(spmdc_sel_pf)])
    dc_output = StructArray(map(VectorOfArrays, columns(dc_output_novv)))

    detdata_output = detdata_output_pf.(detdata)

    return StructVector(merge(columns(cal_output), columns(dc_output), columns(detdata_output)))
end
export calibrate_spm_detector_data


function _single_fiber_esum(
    t_win::AbstractInterval, spm_t::AbstractVector{<:Number},
    spmdc::AbstractVector{<:Number}, spm_pe::AbstractVector{<:Number},
    pe_trig_threshold::Quantity{<:Real}
)
    s::eltype(spm_pe) = zero(eltype(spm_pe))
    for i in eachindex(spm_t)
        if spm_t[i] in t_win && !spmdc[i] && spm_pe[i] > pe_trig_threshold
            s += spm_pe[i]
        end
    end
    return s
end


function _lar_cut(
    colnames::Tuple{Symbol, Symbol, Symbol, Symbol},
    t_win::AbstractInterval, spm_t::AbstractVector{<:AbstractVector{<:Number}},
    spmdc::AbstractVector{<:AbstractVector{<:Number}}, spm_pe::AbstractVector{<:AbstractVector{<:Number}},
    pe_trig_threshold::Quantity{<:Real}, pe_sum_threshold::Quantity{<:Real}, multiplicity_threshold::Int
)
    PEt = eltype(eltype(spm_pe))
    n = length(spm_t)
    pe_sum_per_fiber = Vector{PEt}(undef, n)
    n_over_thresh::Int = 0
    pe_sum::PEt = zero(PEt)
    for i in 1:n
        s_i = _single_fiber_esum(t_win, spm_t[i], spmdc[i], spm_pe[i], pe_trig_threshold)
        pe_sum_per_fiber[i] = s_i
        if s_i > zero(s_i)
            n_over_thresh += 1
        end
        pe_sum += s_i
    end

    lar_cut = Bool(n_over_thresh >= multiplicity_threshold || pe_sum >= pe_sum_threshold)

    return NamedTuple{colnames, Tuple{Bool, PEt, Int, Vector{PEt}}}((lar_cut, pe_sum, n_over_thresh, pe_sum_per_fiber))
end

# Per-fiber window-summed PE (VoV per event). No aggregates and no cut Bool —
# used for the prompt/delayed windows where downstream analysis works on the per-SiPM values directly.
function _window_pe_sum_per_fiber(
    t_win::AbstractInterval, spm_t::AbstractVector{<:AbstractVector{<:Number}},
    spmdc::AbstractVector{<:AbstractVector{<:Number}}, spm_pe::AbstractVector{<:AbstractVector{<:Number}},
    pe_trig_threshold::Quantity{<:Real}
)
    PEt = eltype(eltype(spm_pe))
    n = length(spm_t)
    pe_sum_per_fiber = Vector{PEt}(undef, n)
    for i in 1:n
        pe_sum_per_fiber[i] = _single_fiber_esum(t_win, spm_t[i], spmdc[i], spm_pe[i], pe_trig_threshold)
    end
    return pe_sum_per_fiber
end

# Additional sum-only windows. (config_key, column suffix)
const _EXTRA_LAR_WINDOWS = (
    (:ged_sum_window_prompt,  "prompt"),
    (:ged_sum_window_delayed, "delayed"),
)

function _build_lar_cut(data::LegendData, sel::AnyValiditySelection, global_events::AbstractVector{<:NamedTuple}, geds_t0::AbstractVector{<:Unitful.Time{<:Real}}, e_filter::Symbol)
    dataprod_larcut = get_spms_evt_lar_cut_props(data, sel)
    dataprod_larcut_filter = dataprod_larcut.energy_types[e_filter]

    spm_t =  getproperty(global_events.spms, Symbol(dataprod_larcut_filter.pos))
    spmdc =  getproperty(global_events.spms, Symbol(dataprod_larcut_filter.is_dc))
    spm_pe = getproperty(global_events.spms, e_filter)

    ged_sum_window = dataprod_larcut.ged_sum_window
    t_wins = ClosedInterval.(geds_t0 .+ first(ged_sum_window), geds_t0 .+ last(ged_sum_window))

    pe_trig_threshold = dataprod_larcut_filter.pe_trig_threshold
    pe_sum_threshold = dataprod_larcut.pe_sum_threshold
    multiplicity_threshold = dataprod_larcut.multiplicity_threshold

    colnames = Tuple(Symbol.("$(e_filter)_" .* ["lar_cut", "spms_win_pe_sum", "spms_win_multiplicity", "spms_pe_sum"]))
    return StructArray(_lar_cut.(Ref(colnames), t_wins, spm_t, spmdc, spm_pe, Ref(pe_trig_threshold), Ref(pe_sum_threshold), Ref(multiplicity_threshold)))
end

function _build_window_pe_sum_per_fiber(data::LegendData, sel::AnyValiditySelection, global_events::AbstractVector{<:NamedTuple}, geds_t0::AbstractVector{<:Unitful.Time{<:Real}}, e_filter::Symbol, window_key::Symbol, suffix::String)
    dataprod_larcut = get_spms_evt_lar_cut_props(data, sel)
    dataprod_larcut_filter = dataprod_larcut.energy_types[e_filter]

    spm_t =  getproperty(global_events.spms, Symbol(dataprod_larcut_filter.pos))
    spmdc =  getproperty(global_events.spms, Symbol(dataprod_larcut_filter.is_dc))
    spm_pe = getproperty(global_events.spms, e_filter)

    window = getproperty(dataprod_larcut, window_key)
    t_wins = ClosedInterval.(geds_t0 .+ first(window), geds_t0 .+ last(window))

    pe_trig_threshold = dataprod_larcut_filter.pe_trig_threshold

    colname = Symbol("$(e_filter)_spms_pe_sum_$(suffix)")
    pe_sum_per_fiber = _window_pe_sum_per_fiber.(t_wins, spm_t, spmdc, spm_pe, Ref(pe_trig_threshold))
    return (; colname => pe_sum_per_fiber)
end

function _build_lar_cut(data::LegendData, sel::AnyValiditySelection, global_events::AbstractVector{<:NamedTuple}, geds_t0::AbstractVector{<:Unitful.Time{<:Real}})
    dataprod_larcut = get_spms_evt_lar_cut_props(data, sel)
    energy_types = keys(dataprod_larcut.energy_types)

    is_valid_lar_propfunc = ljl_propfunc(dataprod_larcut.is_valid_lar)

    # LAr cut (full window with cut Bool) — existing behaviour, plus per-fiber pe_sums column
    lar_cut_cols = columns.(_build_lar_cut.(Ref(data), Ref(sel), Ref(global_events), Ref(geds_t0), energy_types))

    # Sum-only windows (prompt, delayed) — only per-fiber pe_sum (VoV)
    extra_cols = NamedTuple[]
    for (wkey, suffix) in _EXTRA_LAR_WINDOWS
        for e_filter in energy_types
            push!(extra_cols, _build_window_pe_sum_per_fiber(data, sel, global_events, geds_t0, e_filter, wkey, suffix))
        end
    end

    lar_cut = StructArray(merge(lar_cut_cols..., extra_cols...))
    result = StructArray(merge((is_valid_lar = is_valid_lar_propfunc.(lar_cut),), columns(lar_cut)))
    # Per-fiber pe_sum columns are Vector{Vector{T}}; convert to VectorOfVectors for proper HDF5 VoV serialization.
    return StructArray(map(_fix_vov, columns(result)))
end
