# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

using LegendEventAnalysis, LegendMakie, CairoMakie, Dates, Unitful, Test
import Makie

@testset "Time-series plotting extension" begin
    @test !isnothing(Base.get_extension(LegendEventAnalysis, :LegendEventAnalysisLegendMakieExt))
    time = DateTime(2026, 1, 1) .+ Day.(0:4)
    pulser = gain_stability(time, [1000, 1001, 999, 1002, 1000]; n_smooth = 3)
    cusp = gain_stability(time, [2000, 2003, 1999, 2001, 2000]; n_smooth = 3)
    fig = lplot(pulser; label = "Pulser", ylabel = "ΔE (keV)", boundaries = time[2:3], figsize = (800, 400),
        watermark = true, watermark_options = (; final = false))
    ax = content(fig[1, 1])
    @test fig isa Makie.Figure
    @test ax.ylabel[] == "ΔE (keV)"
    @test count(p -> p isa Makie.Band, ax.scene.plots) == 2
    @test count(p -> p isa Makie.Lines, ax.scene.plots) == 1
    @test lplot!(ax, cusp; color = :red, label = "E_cusp") === fig
    @test count(p -> p isa Makie.Lines, ax.scene.plots) == 2
    Makie.axislegend(ax)
    mktempdir() do dir
        path = joinpath(dir, "gain.png")
        save(path, fig)
        @test filesize(path) > 0
    end

    sparse = TimeSeriesStruct(time, [1, 2, 1, 3, 2] .* u"keV")
    fig = lplot(sparse; plot = Makie.scatter!, title = "One value per run")
    ax = content(fig[1, 1])
    @test only(ax.scene.plots) isa Makie.Scatter
    @test ax.title[] == "One value per run"
    @test lplot!(pulser; col = 2, sigmas = (), axis = (; limits = (nothing, (-8, 8)))) === fig
    @test only(content(fig[1, 2]).scene.plots) isa Makie.Lines

    unitful = TimeSeriesStruct(time, [1, 2, 1, 3, 2] .* u"keV";
        uncertainty = fill(0.1u"keV", 5))
    fig = lplot(unitful; sigmas = (1,))
    @test count(p -> p isa Makie.Band, content(fig[1, 1]).scene.plots) == 1
end
