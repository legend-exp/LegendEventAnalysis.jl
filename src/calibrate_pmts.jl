# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    calibrate_pmt_channel_data(data::LegendData, sel::ValiditySelection, detector::DetectorId, channel_data::AbstractVector)

Apply the calibration specified by `data` and `sel` for the given PMT
`detector` to the single-channel `channel_data` for that detector.

Also calculates the configured cut/flag values.
"""
function calibrate_pmt_channel_data(data::LegendData, sel::AnyValiditySelection, detector::DetectorId, channel_data::AbstractVector; 
    e_cal_pars_type::Symbol=:rpars, e_cal_pars_cat::Symbol=:pmtcal, 
    keep_chdata::Bool=false)
    chdata = channel_data[:]

    pmtcal_pf = get_pmt_cal_propfunc(data, sel, detector; pars_type=e_cal_pars_type, pars_cat=e_cal_pars_cat)
    cut_is_physical_trig = get_pmt_is_physical_trig_propfunc(data, sel, detector; pars_type=e_cal_pars_type)

    # get additional cols to be parsed into the event tier
    chdata_output_pf = if keep_chdata
        PropSelFunction{propertynames(chdata)}()
    else
        get_pmts_evt_chdata_propfunc(data, sel)
    end

    cal_output = pmtcal_pf.(chdata)

    cal_chdata = StructArray(merge(columns(cal_output), columns(chdata)))

    chdata_output = chdata_output_pf.(chdata)

    is_physical_trig = cut_is_physical_trig.(cal_chdata)
    
    return StructVector(merge(columns(cal_output), columns(chdata_output), columns(is_physical_trig)))
end
export calibrate_pmt_channel_data

function _muon_cut(
    colnames::Tuple{Symbol, Symbol, Symbol},
    is_physical_trig::AbstractVector{<:Bool}, pmt_pe::AbstractVector{<:Number},
    pe_ch_threshold::Quantity{<:Real}, pe_sum_threshold::Quantity{<:Real}, multiplicity_threshold::Int
)   
    trig_ch = findall(is_physical_trig .&& pmt_pe .> pe_ch_threshold)
    pe_sum::eltype(eltype(pmt_pe)) = sum(pmt_pe[trig_ch])
    n_over_thresh::Int = length(trig_ch)

    muon_cut = Bool(n_over_thresh > multiplicity_threshold || pe_sum > pe_sum_threshold)

    return NamedTuple{colnames, Tuple{Bool, eltype(eltype(pmt_pe)), Int}}([muon_cut, pe_sum, n_over_thresh])
end

function _build_muon_cut(data::LegendData, sel::AnyValiditySelection, pmt_events::AbstractVector{<:NamedTuple}, e_filter::Symbol)
    dataprod_muoncut = get_pmts_evt_muon_cut_props(data, sel)

    pmt_pe = getproperty(pmt_events, e_filter)
    is_physical_trig = getproperty(pmt_events, Symbol(replace(string(e_filter), "_cal" => "_is_physical_trig")))

    pe_ch_threshold = dataprod_muoncut.pe_ch_threshold
    pe_sum_threshold = dataprod_muoncut.pe_sum_threshold
    multiplicity_threshold = dataprod_muoncut.multiplicity_threshold

    colnames = Tuple(Symbol.("$(e_filter)_" .* ["muon_cut", "pmts_pe_sum", "pmts_multiplicity"]))
    return StructArray(_muon_cut.(Ref(colnames), is_physical_trig, pmt_pe, Ref(pe_ch_threshold), Ref(pe_sum_threshold), Ref(multiplicity_threshold)))
end

function _build_muon_cut(data::LegendData, sel::AnyValiditySelection, pmt_events::AbstractVector{<:NamedTuple})
    dataprod_muoncut = get_pmts_evt_muon_cut_props(data, sel)
    energy_types = Symbol.(dataprod_muoncut.energy_types)

    is_valid_muon_propfunc = ljl_propfunc(dataprod_muoncut.is_valid_muon)
    
    muon_cut = StructArray(merge(columns.(_build_muon_cut.(Ref(data), Ref(sel), Ref(pmt_events), energy_types))...))
    return StructArray(merge((is_valid_muon = is_valid_muon_propfunc.(muon_cut),), columns(muon_cut)))
end

function _muon_evt_cut(
    t_win::AbstractInterval, pmt_tstart::AbstractVector{<:Number},
    pmt_events_trig::AbstractVector{<:NamedTuple}, empty_evt::NamedTuple
)   
    muon_cut_idx = findall(pmt_tstart .∈ t_win)
    muon_cut = Bool(length(muon_cut_idx) > 0)

    trig_evt::typeof(first(pmt_events_trig)) = if muon_cut
        pmt_events_trig[first(muon_cut_idx)]
    else
        empty_evt
    end

    return merge((is_valid_muon = !muon_cut, n_hit_window = length(muon_cut_idx)), trig_evt)
end


function _build_muon_evt_cut(data::LegendData, sel::AnyValiditySelection, global_events::AbstractVector{<:NamedTuple}, pmt_events::AbstractVector{<:NamedTuple})
    output_colnames = (:timestamp, :pe_cal, :multiplicity)
    if isempty(pmt_events)
        return StructArray(fill(NamedTuple{output_colnames}([0.0u"s", 0.0u"e_au", 0]), length(global_events)))
    elseif count(pmt_events.is_valid_muon) == 0
        return StructArray(fill(NamedTuple{output_colnames}([0.0u"s", 0.0u"e_au", 0]), length(global_events)))
    end
    geds_t0_absolute = global_events.tstart .+ global_events.geds.t0_start
    
    dataprod_muoncut = get_pmts_evt_muon_cut_props(data, sel)
    
    ged_sum_window = dataprod_muoncut.ged_sum_window

    t_wins = ClosedInterval.(geds_t0_absolute .+ first(ged_sum_window), geds_t0_absolute .+ last(ged_sum_window))

    evtdata_output_pf = get_pmts_evt_evtdata_propfunc(data, sel)
    pmt_events_trig = evtdata_output_pf.(pmt_events[findall(pmt_events.is_valid_muon)])

    @assert length(propertynames(pmt_events_trig)) == length(output_colnames) && all(hasproperty.(Ref(pmt_events_trig), output_colnames)) "Invalid `evtdata_output` properties"

    empty_evt::typeof(first(pmt_events_trig)) = typeof(first(pmt_events_trig))([zero(first(pmt_events_trig)[k]) for k in keys(first(pmt_events_trig))])

    return StructArray(_muon_evt_cut.(t_wins, Ref(pmt_events.tstart), Ref(pmt_events_trig), Ref(empty_evt)))
end