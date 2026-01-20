# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).


function _data_with_detector(all_det_data::AbstractDict{<:DetectorIdLike}, det::DetectorIdLike)
    data = all_det_data[det][:]
    det_cols = (detector = fill(DetectorId(det), length(data)), detevtno = collect(eachindex(data)))
    StructVector(merge(det_cols, columns(data)))
end

"""
    function flatten_over_detectors(
        data::AbstractDict{<:DetectorIdLike},
        detectors::AbstractVector{<:DetectorIdLike} = collect(keys(data))
    )

Flatten per-detector data `data` to a single `StructArrays.StructVector` by
concatenating its table-like values and adding the columns `detector` and
`detevtno`.

`data` must a dictionary of in-memory or on-disk table-like objects, keyed by
detector-IDs. It may, e.g. be a `Dict` with values that are
`StructArrays.StructVector`, `TypedTables.Table` or similar, but may also be a
`LegendHDF5IO.LHDataStore`. Note that on-disk data will be read into memory
as a whole.

Returns a `NamedTuple{(:result, :per_detector)}`: `result` is the flattened
data and `per_detector` is a detector-indexed dictionary of views into `result`.
"""
function flatten_over_detectors(
    data::AbstractDict{<:DetectorIdLike},
    detectors::AbstractVector{<:DetectorIdLike} = collect(keys(data))
)
    flat_result = _data_with_detector(data, DetectorId(first(detectors)))
    view_ranges = [firstindex(flat_result):lastindex(flat_result)]
    @showprogress desc="Flattening data over detectors" for det in detectors[begin+1:end]
        chunk = _data_with_detector(data, DetectorId(det))
        push!(view_ranges, lastindex(flat_result)+1:lastindex(flat_result)+length(chunk))
        append!(flat_result, chunk)
    end
    dict_result = Dict(broadcast((det, r) -> DetectorId(det) => view(flat_result, r), detectors, view_ranges))
    return (result = flat_result, per_detector = dict_result)
end
export flatten_over_detectors
