# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

"""
    TimeSeriesStruct(time, values; uncertainty = nothing)

One quantity sampled over time, optionally with one-sigma uncertainties.
Vectors must have matching lengths; timestamps and units are preserved.
Values may be individual samples or averages at their corresponding times.
Load LegendMakie and Makie to plot with `lplot` / `lplot!`.
"""
struct TimeSeriesStruct{T<:AbstractVector,V<:AbstractVector,U<:Union{Nothing,AbstractVector}}
    time::T
    values::V
    uncertainty::U

    function TimeSeriesStruct(time::AbstractVector, values::AbstractVector; uncertainty = nothing)
        length(time) == length(values) || throw(DimensionMismatch("time and values must have equal lengths"))
        isnothing(uncertainty) || length(uncertainty) == length(values) ||
            throw(DimensionMismatch("uncertainty and values must have equal lengths"))
        new{typeof(time),typeof(values),typeof(uncertainty)}(time, values, uncertainty)
    end
end
export TimeSeriesStruct

"""
    gain_stability(time, energy; Qbb = 2039.0, n_ref = 500, n_smooth = 200)

Return a `TimeSeriesStruct` of gain-induced energy shifts at `Qbb` for one
energy estimator. Subtract the median of the first `n_ref` samples and scale
by `Qbb / first(energy)`, then compute a centered rolling mean and sample
standard deviation. Windows are truncated at the edges; a singleton has zero
spread. The window radius is `n_smooth ÷ 2`. Times stay absolute and unchanged.
The caller must provide nonempty, finite inputs, a nonzero first energy, and
positive `n_ref` and `n_smooth`. Remove nonfinite time/energy pairs before calling
this function; it does not filter or validate these inputs.
Use two reports to overlay pulser and reconstructed-energy stability.
"""
function gain_stability(time::AbstractVector, energy::AbstractVector; Qbb = 2039.0, n_ref::Integer = 500, n_smooth::Integer = 200)
    reference = median(@view energy[1:min(n_ref, length(energy))])
    shifts = (energy .- reference) ./ first(energy) .* Qbb
    halfwidth = n_smooth ÷ 2
    windows = [@view shifts[max(1, i-halfwidth):min(end, i+halfwidth)] for i in eachindex(shifts)]
    values = mean.(windows)
    uncertainty = [std(w; corrected = length(w) > 1) for w in windows]
    TimeSeriesStruct(time, values; uncertainty)
end
export gain_stability
