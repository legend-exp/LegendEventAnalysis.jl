# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

using Test
using LegendEventAnalysis
import Documenter

Documenter.DocMeta.setdocmeta!(
    LegendEventAnalysis,
    :DocTestSetup,
    :(using LegendEventAnalysis);
    recursive=true,
)
Documenter.doctest(LegendEventAnalysis)
