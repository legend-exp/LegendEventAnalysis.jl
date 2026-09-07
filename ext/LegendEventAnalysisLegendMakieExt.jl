# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

module LegendEventAnalysisLegendMakieExt

using LegendEventAnalysis: TimeSeriesStruct
import LegendMakie
import Makie

function LegendMakie.lplot!(ax::Makie.Axis, report::TimeSeriesStruct;
    plot = Makie.lines!, color = :blue, label = nothing, sigmas = (3, 1), kwargs...)

    plot(ax, report.time, report.values; color, label, kwargs...)
    if !isnothing(report.uncertainty)
        # Bands use numeric coordinates in the axis's time and value units.
        time = Makie.convert_dim_value.(Ref(ax), 1, report.time)
        for sigma in sigmas
            lower = Makie.convert_dim_value.(Ref(ax), 2, report.values .- sigma .* report.uncertainty)
            upper = Makie.convert_dim_value.(Ref(ax), 2, report.values .+ sigma .* report.uncertainty)
            band = Makie.band!(ax.scene, time, lower, upper;
                color = (color, 0.2 / sigma), dim_conversions = Makie.DimConversions())
            Makie.translate!(band, 0, 0, -1)
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
