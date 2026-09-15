# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

module LegendEventAnalysisLegendMakieExt

using LegendEventAnalysis: TimeEvolution
using IntervalSets: ClosedInterval, leftendpoint, rightendpoint
import LegendMakie
import Makie
import Measurements
using Dates: DateTime

function LegendMakie.lplot!(ax::Makie.Axis, ts::TimeEvolution;
    plot = Makie.lines!, color = :blue, label = nothing, sigmas = (3, 1),
    uncertainty_style = :band, whiskerwidth = 7, kwargs...)

    uncertainty_style in (:band, :bars, :none) || throw(ArgumentError("uncertainty_style must be :band, :bars or :none, got $uncertainty_style"))
    values = Measurements.value.(ts.values)
    errors = Measurements.uncertainty.(ts.values)
    plot(ax, ts.time, values; color, label, kwargs...)
    if uncertainty_style != :none && !all(iszero, errors)
        # Uncertainties use numeric coordinates in the axis's time and value units.
        time = Makie.convert_dim_value.(Ref(ax), 1, ts.time)
        values = Makie.convert_dim_value.(Ref(ax), 2, values)
        errors = Makie.convert_dim_value.(Ref(ax), 2, errors)
        if uncertainty_style == :bars
            Makie.errorbars!(ax.scene, time, values, errors;
                color, whiskerwidth, dim_conversions = Makie.DimConversions())
        else
            for sigma in sigmas
                band = Makie.band!(ax.scene, time, values .- sigma .* errors, values .+ sigma .* errors;
                    color = (color, 0.2 / sigma), dim_conversions = Makie.DimConversions())
                Makie.translate!(band, 0, 0, -1)
            end
        end
    end
    ax.parent
end

function LegendMakie.lplot!(ax::Makie.Axis, ts::TimeEvolution{<:AbstractVector{<:DateTime}, <:AbstractVector{<:ClosedInterval}};
    color = :blue, label = nothing, kwargs...)
    lo, hi = leftendpoint.(ts.values), rightendpoint.(ts.values)
    Makie.lines!(ax, ts.time, lo; color, label, kwargs...)
    Makie.lines!(ax, ts.time, hi; color, kwargs...)
    time = Makie.convert_dim_value.(Ref(ax), 1, ts.time)
    band = Makie.band!(ax.scene, time, Makie.convert_dim_value.(Ref(ax), 2, lo), Makie.convert_dim_value.(Ref(ax), 2, hi);
        color = (color, 0.2), dim_conversions = Makie.DimConversions())
    Makie.translate!(band, 0, 0, -1)
    ax.parent
end

function LegendMakie.lplot!(ts::TimeEvolution;
    row = 1, col = 1, xlabel = "Time", ylabel = "", title = "", ylims = nothing,
    axis = (;), boundaries = [], watermark = false, watermark_options = (;), kwargs...)

    fig = Makie.current_figure()
    ax = Makie.Axis(fig[row, col]; xlabel, ylabel, title, axis...)
    LegendMakie.lplot!(ax, ts; kwargs...)
    isnothing(ylims) || Makie.ylims!(ax, ylims...)
    isempty(boundaries) || Makie.vlines!(ax.scene, Makie.convert_dim_value.(Ref(ax), 1, boundaries);
        color = :gray, linestyle = :dash, dim_conversions = Makie.DimConversions())
    Makie.current_axis!(ax)
    watermark && LegendMakie.add_watermarks!(; watermark_options...)
    fig
end

end
