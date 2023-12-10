# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).


"""
    build_global_event_map(local_events::StructArray; ts_window::Number = 25u"μs")

Build a map of global events based on `local_events`.

`data` must contain columns `channel`, `chevtno` and `timestamp`. It will
typically be the result of [`flatten_over_channels`](@ref).

Returns a `StructArray` that contains the columns `start`, `channels`,
`localevents` and `timestamps`, sorted by `start` globally and along
`timestamps` in each row.

The `start` column contains the start time of each event, the `channels`,
`chevts` and `timestamps` columns are vectors of vectors that contain the
channel-id, per-channel event numbers and per-channel
timestamps that have been associated with each respective events.

Per-channel events are accociated with the same global event if their
timestamps fall within a time windows of length `ts_window`.
"""
function build_global_event_map(data::StructVector; ts_window::Number = 25u"μs")
    tssortidxs = sortperm(data.timestamp)
    local_channel = data.channel[tssortidxs]
    local_chevtno = data.chevtno[tssortidxs]
    local_timestamp = data.timestamp[tssortidxs]

    this_evt_start = zero(eltype(data.timestamp))
    this_evt_dataidxs = Vector{Int}()
    this_evt_channels = Vector{eltype(data.channel)}()
    this_evt_chevtnos = Vector{Int}()
    this_evt_timestamps = Vector{eltype(data.timestamp)}()

    all_start = Vector{eltype(data.timestamp)}()
    all_multiplicity = Vector{Int}()
    all_dataidxs = VectorOfVectors{Int}()
    all_channels = VectorOfVectors{eltype(data.channel)}()
    all_chevtnos = VectorOfVectors{Int}()
    all_timestamps = VectorOfVectors{eltype(data.timestamp)}()

    function _flush_evt(i)
        if !isempty(this_evt_channels)
            push!(all_start, this_evt_start)
            push!(all_dataidxs, this_evt_dataidxs); empty!(this_evt_dataidxs)
            push!(all_multiplicity, length(this_evt_chevtnos))
            push!(all_chevtnos, this_evt_chevtnos); empty!(this_evt_chevtnos)
            push!(all_channels, this_evt_channels); empty!(this_evt_channels)
            push!(all_timestamps, this_evt_timestamps); empty!(this_evt_timestamps)
        end
    end

    first_entry = true
    for i in eachindex(local_timestamp)
        ch, chevtno, ts = local_channel[i], local_chevtno[i], local_timestamp[i]
        if first_entry || ts - this_evt_start > ts_window
            _flush_evt(i)
            this_evt_start = ts
        end
        push!(this_evt_dataidxs, tssortidxs[i])
        push!(this_evt_channels, ch)
        push!(this_evt_chevtnos, chevtno)
        push!(this_evt_timestamps, ts)

        first_entry = false
    end

    _flush_evt(-1)

    return StructVector(
        tstart = all_start, multiplicity = all_multiplicity, dataidx = all_dataidxs,
        chevtno = all_chevtnos, channel = all_channels, timestamp = all_timestamps
    )
end
export build_global_event_map


"""
    apply_event_map(data::StructVector, evtmap::StructVector)

Apply the event map `evtmap` to the data `data`.

`data` will typically be the result of [`flatten_over_channels`](@ref) and
`evtmap` the result of [`build_global_event_map`](@ref).
"""
function apply_event_map(data::StructVector, evtmap::StructVector)
    data_cols = columns(data)
    evt_cols = map(c -> VectorOfVectors{eltype(c)}(), data_cols)

    map(data_cols, evt_cols) do c, ac
        for idx in evtmap.dataidx
            push!(ac, c[idx])
        end
        nothing
    end

    @assert evt_cols.timestamp == evtmap.timestamp

    StructVector(merge(columns(evtmap), evt_cols))
end
export apply_event_map


"""
    function build_global_events(
        data::AbstractDict{<:ChannelIdLike},
        channels::AbstractVector{<:ChannelIdLike} = collect(keys(data));
        ts_window::Number = 25u"μs"
    )

Build global events from a dictionary of per-channel events

Per-channel events are accociated with the same global event if their
timestamps fall within a time windows of length `ts_window`.

`data` must a dictionary of in-memory or on-disk table-like objects, keyed by
channel-IDs. It may, e.g. be a `Dict` with values that are
`StructArrays.StructVector`, `TypedTables.Table` or similar, but may also be a
`LegendHDF5IO.LHDataStore`. Note that on-disk data will be read into memory
as a whole.
"""
function build_global_events(
    data::AbstractDict{<:ChannelIdLike},
    channels::AbstractVector{<:ChannelIdLike} = collect(keys(data));
    ts_window::Number = 25u"μs"
)
    flat_data = flatten_over_channels(data, channels).result
    evtmap = build_global_event_map(flat_data; ts_window = ts_window)
    return apply_event_map(flat_data, evtmap)
end
export build_global_events
