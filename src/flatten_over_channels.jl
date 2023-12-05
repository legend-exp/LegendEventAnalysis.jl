# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).


function _data_with_channel(all_ch_data::AbstractDict{<:ChannelIdLike}, ch::ChannelIdLike)
    data = all_ch_data[ch][:]
    ch_cols = (channel = fill(Int(ch), length(data)), chevtno = collect(eachindex(data)))
    StructVector(merge(ch_cols, columns(data)))
end

"""
    function flatten_over_channels(
        data::AbstractDict{<:ChannelIdLike},
        channels::AbstractVector{<:ChannelIdLike} = collect(keys(data))
    )

Flatten per-channel data `data` to a single `StructArrays.StructVector` by
concatenating its table-like values and adding the columns a `channel` and
`chevtno`.

`data` must a dictionary of in-memory or on-disk table-like objects, keyed by
channel-IDs. It may, e.g. be a `Dict` with values that are
`StructArrays.StructVector`, `TypedTables.Table` or similar, but may also be a
`LegendHDF5IO.LHDataStore`. Note that on-disk data will be read into memory
as a whole.

Returns a `NamedTuple{(:result, :per_channel)}`: `result` is the flattened
data and `per_channel` is a channel-indexed dictionary of views into `result`.
"""
function flatten_over_channels(
    data::AbstractDict{<:ChannelIdLike},
    channels::AbstractVector{<:ChannelIdLike} = collect(keys(data))
)
    flat_result = _data_with_channel(data, ChannelId(first(channels)))
    view_ranges = [firstindex(flat_result):lastindex(flat_result)]
    @showprogress desc="Flattening data over channels" for ch in channels[begin+1:end]
        chunk = _data_with_channel(data, ChannelId(ch))
        push!(view_ranges, lastindex(flat_result)+1:lastindex(flat_result)+length(chunk))
        append!(flat_result, chunk)
    end
    dict_result = Dict(broadcast((ch, r) -> ChannelId(ch) => view(flat_result, r), channels, view_ranges))
    return (result = flat_result, per_channel = dict_result)
end
export flatten_over_channels
