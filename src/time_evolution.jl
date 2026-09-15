# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    TimeEvolution(time, values)

One quantity sampled at absolute times. `time` is stored as `DateTime`;
unitful unix timestamps are converted. `values` is stored as given and may
carry units and/or `Measurement` uncertainties. Timestamps must be sorted and
match the axes of `values`.

Derived series: [`smooth`](@ref), [`relative`](@ref), [`rate`](@ref).
Load LegendMakie and a Makie backend to plot with `lplot` / `lplot!`.
"""
struct TimeEvolution{T<:AbstractVector{<:DateTime}, V<:AbstractVector}
    time::T
    values::V

    function TimeEvolution{T,V}(time, values) where {T,V}
        axes(time) == axes(values) || throw(DimensionMismatch("time and values must share axes, got $(axes(time)) and $(axes(values))"))
        issorted(time) || throw(ArgumentError("timestamps must be sorted"))
        new{T,V}(time, values)
    end
end
export TimeEvolution

function TimeEvolution(time::AbstractVector, values::AbstractVector)
    t = _datetime(time)
    TimeEvolution{typeof(t), typeof(values)}(t, values)
end

Base.length(ts::TimeEvolution) = length(ts.time)

_datetime(t::DateTime) = t
_datetime(t::Quantity{<:Real, Unitful.𝐓}) = unix2datetime(ustrip(u"s", t))
_datetime(time::AbstractVector{<:DateTime}) = time
_datetime(time::AbstractVector) = _datetime.(time)

for (M, f) in ((Statistics, :mean), (Statistics, :std), (Statistics, :median), (Base, :extrema), (Base, :minimum), (Base, :maximum))
    @eval $M.$f(ts::TimeEvolution) = $M.$f(ts.values)
end
for f in (:findmax, :findmin)
    @eval Base.$f(ts::TimeEvolution) = ((v, i) = Base.$f(ts.values); (v, ts.time[i]))
end

"""
    ts[t1 .. t2]

Sub-series with timestamps in the closed interval `t1 .. t2` (`DateTime` or unitful unix times).
"""
function Base.getindex(ts::TimeEvolution, span::ClosedInterval)
    idxs = searchsortedfirst(ts.time, _datetime(leftendpoint(span))):searchsortedlast(ts.time, _datetime(rightendpoint(span)))
    TimeEvolution(ts.time[idxs], ts.values[idxs])
end

"""
    filter(f, ts::TimeEvolution)

Sub-series of the samples whose values satisfy `f`.
"""
Base.filter(f, ts::TimeEvolution) = (keep = findall(f, ts.values); TimeEvolution(ts.time[keep], ts.values[keep]))

# `DateTime` has millisecond resolution, so every window width is a `Millisecond`.
_period(Δt::Dates.Period) = Dates.Millisecond(Δt)
_period(Δt::Quantity{<:Real, Unitful.𝐓}) = Dates.Millisecond(round(Int, ustrip(u"ms", Δt)))

_quantity(Δt::Dates.Period) = Dates.Millisecond(Δt).value * u"ms"
_quantity(Δt::Quantity{<:Real, Unitful.𝐓}) = Δt

# Consecutive half-open windows `[start, start + Δt)` of width `Δt`, starting at
# the first timestamp and covering the last one, as `(start, idxs)` pairs where
# `idxs` are the indices of the timestamps inside the window.
function _windows(time::AbstractVector{<:DateTime}, Δt::Dates.Period)
    Δt > zero(Δt) || throw(ArgumentError("window width must be positive, got $Δt"))
    issorted(time) || throw(ArgumentError("timestamps must be sorted"))
    [(start, searchsortedfirst(time, start):searchsortedfirst(time, start + Δt) - 1)
        for start in first(time):Δt:last(time)]
end

"""
    smooth(ts::TimeEvolution, Δt, f = mean)

Reduce `ts` over consecutive windows of width `Δt` (a `Dates.Period` or unitful
time) with the statistic `f`. Each non-empty window yields one sample at the
window center. `f = mean` gives the window mean with the sample standard
deviation as uncertainty (zero for a single sample; input uncertainties are
not propagated); `f = extrema` gives a `ClosedInterval`; any other `f` is
applied to the window values as is.
"""
function smooth(ts::TimeEvolution, Δt, f = mean)
    Δt = _period(Δt)
    windows = filter(w -> !isempty(last(w)), _windows(ts.time, Δt))
    TimeEvolution([start + Δt ÷ 2 for (start, _) in windows], [f(view(ts.values, idxs)) for (_, idxs) in windows])
end
smooth(ts::TimeEvolution, Δt, ::typeof(mean)) = smooth(ts, Δt, _mean_spread)
smooth(ts::TimeEvolution, Δt, ::typeof(extrema)) = smooth(ts, Δt, v -> ClosedInterval(extrema(v)...))
_mean_spread(v) = (x = Measurements.value.(v); measurement(mean(x), std(x; corrected = length(x) > 1)))
export smooth

"""
    relative(ts::TimeEvolution, reference)

Deviation of `ts` from its initial level in percent. The level is the mean of
the first `reference` samples (an `Integer`) or of the samples within
`reference` (a `Dates.Period` or unitful time) after the first timestamp.
Units cancel; `Measurement` uncertainties propagate.
"""
function relative(ts::TimeEvolution, reference)
    level = mean(view(ts.values, _reference_range(ts, reference)))
    TimeEvolution(ts.time, 100 .* (ts.values .- level) ./ level)
end
export relative

function _reference_range(ts::TimeEvolution, n::Integer)
    0 < n <= length(ts) || throw(ArgumentError("reference sample count must be in 1:$(length(ts)), got $n"))
    firstindex(ts.values):firstindex(ts.values) + n - 1
end
function _reference_range(ts::TimeEvolution, span)
    span = _period(span)
    span > zero(span) || throw(ArgumentError("reference span must be positive, got $span"))
    firstindex(ts.values):searchsortedlast(ts.time, first(ts.time) + span)
end

"""
    rate(timestamps, Δt)

Event rate in Hz over consecutive windows of width `Δt` (a `Dates.Period` or
unitful time), starting at the first timestamp. Only windows fully covered by
the data span are reported, each at its center with a Poisson (`√N`)
uncertainty. `timestamps` are `DateTime` or unitful unix times and must be sorted.
"""
function rate(timestamps::AbstractVector, Δt)
    time = _datetime(timestamps)
    period = _period(Δt)
    windows = filter(w -> first(w) + period <= last(time), _windows(time, period))
    isempty(windows) && throw(ArgumentError("no window of width $Δt fits into the data span $(first(time)) – $(last(time))"))
    counts = [length(idxs) for (_, idxs) in windows]
    TimeEvolution(
        [start + period ÷ 2 for (start, _) in windows],
        uconvert.(u"Hz", measurement.(counts, sqrt.(counts)) ./ _quantity(Δt))
    )
end
export rate
