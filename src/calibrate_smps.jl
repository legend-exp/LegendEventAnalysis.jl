# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_spm_detector_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, detector_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for the given SiPM
`detector` to the single-detector `detector_data` for that detector.

Also calculates the configured cut/flag values and drops the triggers below
the noise threshold of the detector.
"""
function calibrate_spm_detector_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, detector_data::AbstractVector;
    e_cal_pars_type::Symbol=:ppars, e_cal_pars_cat::Symbol=:sipmcal, dc_cut_pars_type::Symbol=:ppars,
    keep_detdata::Bool=false)
    detdata = detector_data[:]

    spmcal_pf = get_spm_cal_propfunc(data, sel, detector; pars_type=e_cal_pars_type, pars_cat=e_cal_pars_cat)
    spmdc_sel_pf = get_spm_dc_sel_propfunc(data, sel, detector; pars_type=dc_cut_pars_type)
    spmdc_cal_pf = get_spm_dc_cal_propfunc(data, sel, detector; pars_type=dc_cut_pars_type)
    spmnoise_pf = get_spm_noise_cut_propfunc(data, sel, detector; pars_type=e_cal_pars_type, pars_cat=e_cal_pars_cat)

    # get additional cols to be parsed into the event tier
    detdata_output_pf = if keep_detdata
        PropSelFunction(propertynames(detdata)...)
    else
        get_spms_evt_detdata_propfunc(data, sel)
    end

    cal_output_novv = spmcal_pf.(detdata)
    cal_output = StructArray(map(VectorOfArrays, columns(cal_output_novv)))

    dc_output_novv = NamedTuple{keys(spmdc_sel_pf)}([spmdc_cal_pf[e_type].(spmdc_sel_pf[e_type].(detdata)) for e_type in keys(spmdc_sel_pf)])
    dc_output = StructArray(map(VectorOfArrays, columns(dc_output_novv)))

    detdata_output = detdata_output_pf.(detdata)

    output = StructVector(merge(columns(cal_output), columns(dc_output), columns(detdata_output)))

    noise_output_novv = spmnoise_pf.(output)
    noise_output = StructArray(map(VectorOfArrays, columns(noise_output_novv)))

    return StructVector(merge(columns(output), columns(noise_output)))
end
export calibrate_spm_detector_data


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


function _spms_win_sum(
    colnames::Tuple{Symbol, Symbol},
    t_win::AbstractInterval, spm_t::AbstractVector{<:AbstractVector{<:Number}},
    spmdc::AbstractVector{<:AbstractVector{<:Number}}, spm_pe::AbstractVector{<:AbstractVector{<:Number}}
)
    n_over_thresh::Int = 0
    pe_sum::eltype(eltype(spm_pe)) = zero(eltype(eltype(spm_pe)))
    for i in eachindex(spm_t)
        s_i = _single_fiber_esum(t_win, spm_t[i], spmdc[i], spm_pe[i])
        if s_i > zero(s_i)
            n_over_thresh += 1
        end
        pe_sum += s_i
    end

    return NamedTuple{colnames, Tuple{eltype(eltype(spm_pe)), Int}}([pe_sum, n_over_thresh])
end

_spm_win_pe_sums(t_win::AbstractInterval, spm_t::AbstractVector, spmdc::AbstractVector, spm_pe::AbstractVector) = [_single_fiber_esum(t_win, spm_t[i], spmdc[i], spm_pe[i]) for i in eachindex(spm_t)]

function _spms_win_sums(windows, t_wins, prefix::String, spm_t::AbstractVector, spmdc::AbstractVector, spm_pe::AbstractVector)
    cols = NamedTuple()
    for win_name in windows
        colnames = (Symbol("$(prefix)$(win_name)_pe_sum"), Symbol("$(prefix)$(win_name)_multiplicity"))
        cols = merge(cols, columns(StructArray(_spms_win_sum.(Ref(colnames), t_wins(win_name), spm_t, spmdc, spm_pe))))
    end
    return cols
end


function _lar_cut(
    t_win::AbstractInterval, spm_t::AbstractVector{<:AbstractVector{<:Number}},
    spmdc::AbstractVector{<:AbstractVector{<:Number}}, spm_pe::AbstractVector{<:AbstractVector{<:Number}},
    pe_sum_threshold::Quantity{<:Real}, multiplicity_threshold::Int
)
    pe_sum, n_over_thresh = _spms_win_sum((:pe_sum, :multiplicity), t_win, spm_t, spmdc, spm_pe)
    return Bool(n_over_thresh >= multiplicity_threshold || pe_sum >= pe_sum_threshold)
end


function _geom_spms(data::LegendData, sel::AnyValiditySelection, dphi_max)
    geometry = data.metadata.hardware.geometry
    spms_geometry = geometry.spms(sel)
    return Dict(
        DetectorId(hpge) => Set(DetectorId(spm) for (spm, rel) in rel_spms if abs(rel.dphi) <= dphi_max[Symbol(spms_geometry[spm].barrel)])
        for (hpge, rel_spms) in geometry.relative(sel)
    )
end

_geom_spms_idxs(geom_spms::AbstractDict, trig_e_det::AbstractVector{<:DetectorId}, spm_det::AbstractVector{<:DetectorId}) = findall(spm -> any(spm in geom_spms[det] for det in trig_e_det), spm_det)


function _build_lar_cut(
    data::LegendData, sel::AnyValiditySelection, global_events::AbstractVector{<:NamedTuple},
    geds_t0::AbstractVector{<:Unitful.Time{<:Real}}, trig_e_det::Union{Nothing, AbstractVector{<:AbstractVector{<:DetectorId}}}, e_filter::Symbol
)
    dataprod_larcut = get_spms_evt_lar_cut_props(data, sel)
    dataprod_larcut_filter = dataprod_larcut.energy_types[e_filter]

    spm_t =  getproperty(global_events.spms, Symbol(dataprod_larcut_filter.pos))
    spmdc =  getproperty(global_events.spms, Symbol(dataprod_larcut_filter.is_dc))
    spm_pe = getproperty(global_events.spms, e_filter)

    windows = sort!(collect(keys(dataprod_larcut.windows)))
    t_wins(win_name) = let sum_window = dataprod_larcut.windows[win_name]
        ClosedInterval.(geds_t0 .+ first(sum_window), geds_t0 .+ last(sum_window))
    end

    # pe sums and multiplicities over all SiPMs for every window
    cols = _spms_win_sums(windows, t_wins, "$(e_filter)_spms_", spm_t, spmdc, spm_pe)

    # cut flags, a cut with dphi_max runs on the SiPMs near the triggered HPGe detectors and stores its own sums
    for cut_name in sort!(collect(keys(dataprod_larcut.cuts)))
        cut = dataprod_larcut.cuts[cut_name]
        if haskey(cut, :dphi_max)
            isnothing(trig_e_det) && continue
            geom_idxs = _geom_spms_idxs.(Ref(_geom_spms(data, sel, cut.dphi_max)), trig_e_det, global_events.spms.detector)
            cut_t, cut_dc, cut_pe = view.(spm_t, geom_idxs), view.(spmdc, geom_idxs), view.(spm_pe, geom_idxs)
            cols = merge(cols, _spms_win_sums(windows, t_wins, "$(e_filter)_spms_$(cut_name)_", cut_t, cut_dc, cut_pe))
        else
            cut_t, cut_dc, cut_pe = spm_t, spmdc, spm_pe
        end
        lar_cut = _lar_cut.(t_wins(Symbol(cut.window)), cut_t, cut_dc, cut_pe, Ref(cut.pe_sum_threshold), Ref(cut.multiplicity_threshold))
        cols = merge(cols, NamedTuple{(Symbol("$(e_filter)_$(cut_name)_lar_cut"),)}((lar_cut,)))
    end

    # pe sums per SiPM for the configured windows, aligned with spms.detector
    spm_cols = NamedTuple()
    for win_name in Symbol.(dataprod_larcut.spm_windows)
        spm_cols = merge(spm_cols, NamedTuple{(Symbol("$(e_filter)_$(win_name)_pe_sum"),)}((VectorOfVectors(_spm_win_pe_sums.(t_wins(win_name), spm_t, spmdc, spm_pe)),)))
    end
    return (cols = cols, spm_cols = spm_cols)
end

function _build_lar_cut(
    data::LegendData, sel::AnyValiditySelection, global_events::AbstractVector{<:NamedTuple},
    geds_t0::AbstractVector{<:Unitful.Time{<:Real}}, trig_e_det::Union{Nothing, AbstractVector{<:AbstractVector{<:DetectorId}}}
)
    dataprod_larcut = get_spms_evt_lar_cut_props(data, sel)
    energy_types = keys(dataprod_larcut.energy_types)

    is_valid_lar_propfunc = ljl_propfunc(dataprod_larcut.is_valid_lar)
    
    lar_cuts = _build_lar_cut.(Ref(data), Ref(sel), Ref(global_events), Ref(geds_t0), Ref(trig_e_det), energy_types)
    lar_cut = StructArray(merge(getproperty.(lar_cuts, :cols)...))
    spm = StructArray(merge((detector = global_events.spms.detector,), getproperty.(lar_cuts, :spm_cols)...))
    return StructArray(merge((is_valid_lar = is_valid_lar_propfunc.(lar_cut),), columns(lar_cut), (spm = spm,)))
end
