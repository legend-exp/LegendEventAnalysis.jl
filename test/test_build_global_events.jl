# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

using LegendEventAnalysis
using Test

using Unitful

@testset "build_global_events" begin
    @testset "flag_coincidences" begin
        timestamps = [0.9, 1.9, 3.5, 4.1, 4.8, 5.1, 6.3, 6.9, 7.6, 8.1] * u"μs"
        ref_timestamps = [2.0, 4.0, 5.2, 7.0] * u"μs"
        expected_flags = [false, true, false, true, false, true, false, true, false, false]
        @test @inferred(flag_coincidences(timestamps, ref_timestamps, ts_window = 0.11u"μs")) == expected_flags
    end
end
