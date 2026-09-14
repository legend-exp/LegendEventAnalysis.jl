# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

module LegendEventAnalysisLegendMakieExt

using LegendEventAnalysis: TimeSeriesStruct
import LegendMakie
import Makie

function LegendMakie.lplot!(ax::Makie.Axis, report::TimeSeriesStruct;
    plot = Makie.lines!, color = :blue, label = nothing, sigmas = (3, 1),
    uncertainty_style = :band, whiskerwidth = 7, kwargs...)

    uncertainty_style in (:band, :bars, :none) || throw(ArgumentError("uncertainty_style must be :band, :bars or :none"))
    plot(ax, report.time, report.values; color, label, kwargs...)
    if !isnothing(report.uncertainty) && uncertainty_style != :none
        # Uncertainties use numeric coordinates in the axis's time and value units.
        time = Makie.convert_dim_value.(Ref(ax), 1, report.time)
        values = Makie.convert_dim_value.(Ref(ax), 2, report.values)
        errors = Makie.convert_dim_value.(Ref(ax), 2, report.uncertainty)
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

function LegendMakie.lplot!(report::TimeSeriesStruct;
    row = 1, col = 1, xlabel = "Time", ylabel = "", title = "",
    axis = (;), boundaries = [], watermark = false, watermark_options = (;), kwargs...)

    fig = Makie.current_figure()
    ax = Makie.Axis(fig[row, col]; xlabel, ylabel, title, axis...)
    LegendMakie.lplot!(ax, report; kwargs...)
    isempty(boundaries) || Makie.vlines!(ax.scene, Makie.convert_dim_value.(Ref(ax), 1, boundaries);
        color = :gray, linestyle = :dash, dim_conversions = Makie.DimConversions())
    Makie.current_axis!(ax)
    watermark && LegendMakie.add_watermarks!(; watermark_options...)
    fig
end

end
