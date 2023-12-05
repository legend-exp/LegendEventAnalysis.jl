# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

__precompile__(true)

module LegendEventAnalysis

using ArraysOfArrays
using StructArrays
using Unitful

using ProgressMeter: @showprogress
using PropertyFunctions: PropertyFunction, @pf

using LegendDataTypes: fast_flatten

include("build_global_events.jl")

end # module
