# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import LegendEventAnalysis

Test.@testset "Aqua tests" begin
    Aqua.test_all(LegendEventAnalysis)
end # testset
