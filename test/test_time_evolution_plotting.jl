# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

using LegendEventAnalysis, LegendMakie, CairoMakie, Dates, Measurements, Unitful, Test
import Makie

@testset "Time-series plotting extension" begin
    @test !isnothing(Base.get_extension(LegendEventAnalysis, :LegendEventAnalysisLegendMakieExt))
    t0 = DateTime(2026, 1, 1)
    time = t0 .+ Minute.(0:99)
    pulser = relative(smooth(TimeEvolution(time, (1000.0 .+ sin.(0:99)) .* u"keV"), 10u"minute"), 20u"minute")
    cusp = relative(smooth(TimeEvolution(time, (2000.0 .+ cos.(0:99)) .* u"keV"), 10u"minute"), 20u"minute")
    fig = lplot(pulser; label = "Pulser", ylabel = "ΔE (%)", boundaries = time[30:40:70], figsize = (800, 400),
        watermark = true, watermark_options = (; final = false))
    ax = content(fig[1, 1])
    @test fig isa Makie.Figure
    @test ax.ylabel[] == "ΔE (%)"
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

    # Values without uncertainties draw no bands.
    sparse = TimeEvolution(time[1:5], [1, 2, 1, 3, 2] .* u"keV")
    fig = lplot(sparse; plot = Makie.scatter!, title = "One value per run")
    ax = content(fig[1, 1])
    @test only(ax.scene.plots) isa Makie.Scatter
    @test ax.title[] == "One value per run"
    @test lplot!(pulser; col = 2, sigmas = (), ylims = (-8, 8)) === fig
    @test only(content(fig[1, 2]).scene.plots) isa Makie.Lines
    @test content(fig[1, 2]).limits[][2] == (-8, 8)

    unitful = TimeEvolution(time[1:5], ([1, 2, 1, 3, 2] .± 0.1) .* u"keV")
    fig = lplot(unitful; sigmas = (1,), ylims = (0u"keV", 4u"keV"))
    @test count(p -> p isa Makie.Band, content(fig[1, 1]).scene.plots) == 1
    @test content(fig[1, 1]).limits[][2] == (0, 4)

    fig = lplot(unitful; plot = Makie.scatterlines!, uncertainty_style = :bars, linestyle = :dash)
    plots = content(fig[1, 1]).scene.plots
    @test count(p -> p isa Makie.Errorbars, plots) == 1
    @test count(p -> p isa Makie.ScatterLines, plots) == 1
    @test !any(p -> p isa Makie.Band, plots)

    fig = lplot(unitful; plot = Makie.scatter!, uncertainty_style = :none)
    @test only(content(fig[1, 1]).scene.plots) isa Makie.Scatter
    @test_throws ArgumentError lplot(unitful; uncertainty_style = :invalid)

    # Interval-valued series (e.g. from `smooth(ts, Δt, extrema)`) draw an envelope.
    env = smooth(TimeEvolution(time, (1000.0 .+ sin.(0:99)) .* u"keV"), 10u"minute", extrema)
    fig = lplot(env; label = "min–max")
    plots = content(fig[1, 1]).scene.plots
    @test count(p -> p isa Makie.Band, plots) == 1 && count(p -> p isa Makie.Lines, plots) == 2

    # Event rates carry Hz units on the value axis.
    r = rate(t0 .+ Second.(0:7199), 10u"minute")
    fig = lplot(r; plot = Makie.scatter!, uncertainty_style = :bars, ylabel = "Rate")
    @test count(p -> p isa Makie.Errorbars, content(fig[1, 1]).scene.plots) == 1
end
