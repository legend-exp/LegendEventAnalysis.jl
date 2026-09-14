# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

using LegendEventAnalysis, Dates, Statistics, Unitful, Test

@testset "Time series" begin
    time = DateTime(2026, 1, 1) .+ Second.(0:4)
    energy = [1000.0, 1002.0, 1004.0, 1006.0, 1008.0]
    report = TimeSeriesStruct(time, energy .* u"keV")
    @test report.time === time
    @test report.values == energy .* u"keV"
    @test isnothing(report.uncertainty)
    @test_throws DimensionMismatch TimeSeriesStruct(time, energy[1:2])
    @test_throws DimensionMismatch TimeSeriesStruct(time, energy; uncertainty = [1.0])

    gain = gain_stability(time, energy; Qbb = 2000.0, n_ref = 2, n_smooth = 3)
    @test gain.time === time
    @test gain.values ≈ [0, 2, 6, 10, 12]
    @test gain.uncertainty ≈ [sqrt(8), 4, 4, 4, sqrt(8)]

    # Preserve the original gain recipe's normalization, including short runs.
    pct = (energy .- first(energy)) ./ first(energy) .* 100
    pct .-= median(pct)
    short = gain_stability(time, energy; n_ref = 500, n_smooth = 1)
    @test short.values ≈ pct ./ 100 .* 2039
    @test all(iszero, short.uncertainty)
    single = gain_stability(time[1:1], energy[1:1])
    @test single.values == single.uncertainty == [0.0]
    units = gain_stability(time, energy .* u"keV"; Qbb = 2000u"keV", n_ref = 2, n_smooth = 3)
    @test units.values ≈ gain.values .* u"keV"
    @test units.uncertainty ≈ gain.uncertainty .* u"keV"
    @test gain_stability(time, energy; Qbb = 2000.0, n_ref = 2, n_smooth = 2).values == gain.values
    @test_throws DimensionMismatch gain_stability(time[1:2], energy)
    # The caller removes nonfinite pairs before computing gain stability.
    raw_time = [0.0, NaN, 2.0, 3.0, Inf]
    raw_energy = [1000.0, 1001.0, NaN, 1002.0, 1003.0]
    finite = isfinite.(raw_time) .& isfinite.(raw_energy)
    filtered = gain_stability(raw_time[finite], raw_energy[finite];
        Qbb = 2000.0, n_ref = 1, n_smooth = 1)
    @test filtered.time == [0.0, 3.0]
    @test filtered.values ≈ [0.0, 4.0]
end
