# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

import Test

Test.@testset "Package LegendEventAnalysis" begin
    # include("test_aqua.jl")

    include("test_build_global_events.jl")
    include("test_lar_cut.jl")
    include("test_time_evolution.jl")
    # Makie fails to build its glyph texture atlas on 32-bit systems
    if Sys.WORD_SIZE == 64
        include("test_time_evolution_plotting.jl")
    end

    include("test_docs.jl")
    isempty(Test.detect_ambiguities(LegendEventAnalysis))
end # testset
