# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_spm_channel_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, channel_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for the given SiPM
`detector` to the single-channel `channel_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_spm_channel_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, channel_data::AbstractVector;
    e_cal_pars_type::Symbol=:ppars, e_cal_pars_cat::Symbol=:sipmcal, dc_cut_pars_type::Symbol=:ppars,
    keep_chdata::Bool=false)
    chdata = channel_data[:]

    spmcal_pf = get_spm_cal_propfunc(data, sel, detector; pars_type=e_cal_pars_type, pars_cat=e_cal_pars_cat)
    spmdc_sel_pf = get_spm_dc_sel_propfunc(data, sel, detector; pars_type=dc_cut_pars_type)
    spmdc_cal_pf = get_spm_dc_cal_propfunc(data, sel, detector; pars_type=dc_cut_pars_type)

    # get additional cols to be parsed into the event tier
    chdata_output_pf = if keep_chdata
        PropSelFunction{propertynames(chdata)}()
    else
        get_spms_evt_chdata_propfunc(data, sel)
    end

    cal_output_novv = spmcal_pf.(chdata)
    cal_output = StructArray(map(VectorOfArrays, columns(cal_output_novv)))

    dc_output_novv = NamedTuple{keys(spmdc_sel_pf)}([spmdc_cal_pf[e_type].(spmdc_sel_pf[e_type].(chdata)) for e_type in keys(spmdc_sel_pf)])
    dc_output = StructArray(map(VectorOfArrays, columns(dc_output_novv)))

    chdata_output = chdata_output_pf.(chdata)

    return StructVector(merge(columns(cal_output), columns(dc_output), columns(chdata_output)))
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


function _lar_cut(
    colnames::Tuple{Symbol, Symbol, Symbol},
    t_win::AbstractInterval, spm_t::AbstractVector{<:AbstractVector{<:Number}},
    spmdc::AbstractVector{<:AbstractVector{<:Number}}, spm_pe::AbstractVector{<:AbstractVector{<:Number}},
    pe_ch_threshold::Quantity{<:Real}, pe_sum_threshold::Quantity{<:Real}, multiplicity_threshold::Int
)
    n_over_thresh::Int = 0
    pe_sum::eltype(eltype(spm_pe)) = zero(eltype(eltype(spm_pe)))
    for i in eachindex(spm_t)
        s_i = _single_fiber_esum(t_win, spm_t[i], spmdc[i], spm_pe[i])
        if s_i > pe_ch_threshold
            n_over_thresh += 1
        end
        pe_sum += s_i
    end

    lar_cut = Bool(n_over_thresh >= multiplicity_threshold || pe_sum >= pe_sum_threshold)

    return NamedTuple{colnames, Tuple{Bool, eltype(eltype(spm_pe)), Int}}([lar_cut, pe_sum, n_over_thresh])
end

function _build_lar_cut(data::LegendData, sel::AnyValiditySelection, global_events::AbstractVector{<:NamedTuple}, e_filter::Symbol)
    geds_t0 = global_events.geds.t0_start

    dataprod_larcut = get_spms_evt_lar_cut_props(data, sel)
    dataprod_larcut_filter = dataprod_larcut.energy_types[e_filter]

    spm_t =  getproperty(global_events.spms, Symbol(dataprod_larcut_filter.pos))
    spmdc =  getproperty(global_events.spms, Symbol(dataprod_larcut_filter.is_dc))
    spm_pe = getproperty(global_events.spms, e_filter)

    ged_sum_window = dataprod_larcut.ged_sum_window

    t_wins = ClosedInterval.(geds_t0 .+ first(ged_sum_window), geds_t0 .+ last(ged_sum_window))

    pe_ch_threshold = dataprod_larcut.pe_ch_threshold
    pe_sum_threshold = dataprod_larcut.pe_sum_threshold
    multiplicity_threshold = dataprod_larcut.multiplicity_threshold

    colnames = Tuple(Symbol.("$(e_filter)_" .* ["lar_cut", "spms_win_pe_sum", "spms_win_multiplicity"]))
    return StructArray(_lar_cut.(Ref(colnames), t_wins, spm_t, spmdc, spm_pe, Ref(pe_ch_threshold), Ref(pe_sum_threshold), Ref(multiplicity_threshold)))
end

function _build_lar_cut(data::LegendData, sel::AnyValiditySelection, global_events::AbstractVector{<:NamedTuple})
    dataprod_larcut = get_spms_evt_lar_cut_props(data, sel)
    energy_types = keys(dataprod_larcut.energy_types)

    is_valid_lar_propfunc = ljl_propfunc(dataprod_larcut.is_valid_lar)
    
    lar_cut = StructArray(merge(columns.(_build_lar_cut.(Ref(data), Ref(sel), Ref(global_events), energy_types))...))
    return StructArray(merge((is_valid_lar = is_valid_lar_propfunc.(lar_cut),), columns(lar_cut)))
end
