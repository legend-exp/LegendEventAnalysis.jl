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
    all_dataidxs = VectorOfVectors{Int}()
    all_channels = VectorOfVectors{eltype(data.channel)}()
    all_chevtnos = VectorOfVectors{Int}()
    all_timestamps = VectorOfVectors{eltype(data.timestamp)}()

    function _flush_evt(i)
        if !isempty(this_evt_channels)
            push!(all_start, this_evt_start)
            push!(all_dataidxs, this_evt_dataidxs); empty!(this_evt_dataidxs)
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
        tstart = all_start, dataidx = all_dataidxs,
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
            push!(ac, view(c, idx))
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


"""
    function build_cross_system_events(
        data::NamedTuple,
        ts_window::Number = 25u"μs"
    )

Build cross-system events.

`data` must be a NamedTuple with properties that represent the names of
experiment (sub)-systems and values that are the result of
[`build_global_events`](@ref) for each system.

Note: Currently requires the `tstart` columns of all systems to be identical.
"""
function build_cross_system_events(data::NamedTuple)
    tstarts = map(x -> x.tstart, data)
    tstart_ref = first(tstarts)
    for i in eachindex(tstarts)[2:end]
        if (tstarts[i] != tstart_ref) throw(ArgumentError("tstart doesn't match across subsystems")) end
    end
    new_cols = merge((tstart = tstart_ref,), data)
    return StructVector(new_cols)
end
export build_cross_system_events

#_similar_empty(x::Number) = typeof(x)(NaN)
#_similar_empty(x::Integer) = typeof(x)(0)
#_similar_empty(V::AbstractVector) = similar(V, 0)


"""
    flag_coincidences(
        timestamps::AbstractVector{<:RealQuantity}, ref_timestamps::AbstractVector{<:RealQuantity};
        ts_window::Number = 125u"μs"
    )

Flag coincidences in `timestamps` with respect to `ref_timestamps`.

Return a boolean vector of the same length as `timestamps` that is `true`
where a timestamp is within `ts_window` of an element of `ref_timestamps` and
`false` otherwise.
"""
function flag_coincidences(
    timestamps::AbstractVector{<:RealQuantity}, ref_timestamps::AbstractVector{<:RealQuantity};
    ts_window::Number = 125u"μs"
)
    flags = similar(timestamps, Bool)
    @assert axes(flags) == axes(timestamps)
    j = firstindex(ref_timestamps)
    last_j = lastindex(ref_timestamps)
    ref_ts_a = ref_timestamps[j]
    j = j < last_j ? j += 1 : j
    ref_ts_b = ref_timestamps[j]
    for i in eachindex(timestamps)
        ts = timestamps[i]
        delta_t_a = abs(ts - ref_ts_a)
        delta_t_b = abs(ts - ref_ts_b)
        if delta_t_b < delta_t_a
            ref_ts_a = ref_ts_b
            j = j < last_j ? j += 1 : j
            ref_ts_b = ref_timestamps[j]
        end
        flags[i] = abs(ts - ref_ts_a) < ts_window
    end
    return flags
end
export flag_coincidences

