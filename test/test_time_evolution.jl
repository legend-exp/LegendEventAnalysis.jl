# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

using LegendEventAnalysis, Dates, Measurements, Statistics, Unitful, Test
using IntervalSets: (..)

@testset "Time series" begin
    t0 = DateTime(2026, 1, 1)
    time = t0 .+ Minute.(0:9)
    energy = collect(1000.0:2.0:1018.0) .* u"keV"

    @testset "Construction" begin
        ts = TimeEvolution(time, energy)
        @test ts.time === time
        @test ts.values === energy
        @test length(ts) == 10

        # Unitful unix timestamps are converted to DateTime.
        unix = datetime2unix.(time) .* u"s"
        @test TimeEvolution(unix, energy).time == time
        @test TimeEvolution(uconvert.(u"ms", unix), energy).time == time

        @test_throws DimensionMismatch TimeEvolution(time, energy[1:2])
        @test_throws "sorted" TimeEvolution(reverse(time), energy)
        @test_throws MethodError TimeEvolution(datetime2unix.(time), energy)
    end

    @testset "smooth" begin
        ts = TimeEvolution(time, energy)
        s = smooth(ts, 4u"minute")
        # Windows [0,4), [4,8), [8,12) minutes; centers at 2, 6, 10 minutes.
        @test s.time == t0 .+ Minute.([2, 6, 10])
        @test Measurements.value.(s.values) ≈ [1003.0, 1011.0, 1017.0] .* u"keV"
        @test Measurements.uncertainty.(s.values) ≈ [std(1000.0:2.0:1006.0), std(1008.0:2.0:1014.0), std(1016.0:2.0:1018.0)] .* u"keV"
        @test smooth(ts, Minute(4)).values == s.values

        by_count = smooth(ts, 4)
        @test by_count.time == t0 .+ Millisecond.([90_000, 330_000, 510_000])
        @test Measurements.value.(by_count.values) == Measurements.value.(s.values)
        @test Measurements.uncertainty.(by_count.values) == Measurements.uncertainty.(s.values)

        @test Measurements.uncertainty(only(smooth(ts, 20u"minute").values)) ≈ std(ustrip.(energy)) * u"keV"
        # A window holding a single sample has zero spread.
        @test all(iszero, Measurements.uncertainty.(smooth(ts, 1u"minute").values))

        # Input uncertainties are not propagated; the spread of the values is used.
        withunc = TimeEvolution(time, energy .± 0.5u"keV")
        @test smooth(withunc, 4u"minute").values == s.values

        @test_throws "positive" smooth(ts, 0u"minute")
        @test_throws "positive" smooth(ts, -1u"minute")

        # Any statistic can be applied per window.
        @test smooth(ts, 4u"minute", maximum).values == [1006.0, 1014.0, 1018.0] .* u"keV"
        @test smooth(ts, 4u"minute", median).values == [1003.0, 1011.0, 1017.0] .* u"keV"
        @test smooth(ts, 4u"minute", std).values ≈ Measurements.uncertainty.(s.values)
        @test smooth(ts, 4u"minute", extrema).values == [1000.0u"keV" .. 1006.0u"keV", 1008.0u"keV" .. 1014.0u"keV", 1016.0u"keV" .. 1018.0u"keV"]
    end

    @testset "Reductions and selection" begin
        ts = TimeEvolution(time, energy)
        @test mean(ts) == 1009.0u"keV"
        @test std(ts) == std(energy)
        @test median(ts) == 1009.0u"keV"
        @test extrema(ts) == (1000.0u"keV", 1018.0u"keV")
        @test (minimum(ts), maximum(ts)) == extrema(ts)
        @test findmax(ts) == (1018.0u"keV", time[end])
        @test findmin(ts) == (1000.0u"keV", time[1])

        sub = ts[time[3] .. time[5]]
        @test sub.time == time[3:5] && sub.values == energy[3:5]
        @test ts[(t0 + Minute(2)) .. (t0 + Minute(4) + Second(30))].time == time[3:5]
        @test ts[datetime2unix(time[3])u"s" .. datetime2unix(time[5])u"s"].time == time[3:5]
        @test isempty(ts[(t0 + Second(1)) .. (t0 + Second(2))].time)

        out = filter(>(1010.0u"keV"), ts)
        @test out.time == time[7:end] && out.values == energy[7:end]
    end

    @testset "relative" begin
        ts = TimeEvolution(time, energy)
        r = relative(ts, 2)
        @test r.time === time
        @test r.values ≈ 100 .* (ustrip.(energy) .- 1001.0) ./ 1001.0
        @test eltype(r.values) <: Real

        # A time span selects the samples within [t0, t0 + span].
        @test relative(ts, 1u"minute").values ≈ r.values
        @test relative(ts, Minute(1)).values ≈ r.values
        @test relative(ts, 30u"s").values ≈ relative(ts, 1).values

        # Uncertainties propagate through the normalization.
        withunc = TimeEvolution(time, energy .± 1.0u"keV")
        rel_unc = relative(withunc, 1)
        @test Measurements.value.(rel_unc.values) ≈ relative(ts, 1).values
        @test eltype(rel_unc.values) <: Measurement

        @test_throws ArgumentError relative(ts, 0)
        @test_throws ArgumentError relative(ts, 11)
        @test_throws "positive" relative(ts, 0u"s")
    end

    @testset "rate" begin
        # 3 events/minute for 10 minutes, then a lone event closing the span.
        events = t0 .+ Second.(0:20:599)
        push!(events, t0 + Minute(10))
        r = rate(events, 2u"minute")
        @test r.time == t0 .+ Minute.([1, 3, 5, 7, 9])
        @test Measurements.value.(r.values) ≈ fill(0.05, 5) .* u"Hz"
        @test Measurements.uncertainty.(r.values) ≈ fill(sqrt(6) / 120, 5) .* u"Hz"
        @test rate(datetime2unix.(events) .* u"s", Minute(2)).values == r.values

        # Only windows fully covered by the data span are reported.
        @test length(rate(events, 3u"minute").time) == 3
        @test_throws "fits" rate(events, 11u"minute")
        @test_throws "sorted" rate(reverse(events), 1u"minute")
    end
end
