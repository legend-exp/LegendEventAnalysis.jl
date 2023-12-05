# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

import Test

Test.@testset "Package LegendEventAnalysis" begin
    # include("test_aqua.jl")

    include("test_build_global_events.jl")

    include("test_docs.jl")
    isempty(Test.detect_ambiguities(LegendEventAnalysis))
end # testset
