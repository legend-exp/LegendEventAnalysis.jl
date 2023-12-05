# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

__precompile__(true)

module LegendEventAnalysis

using ArraysOfArrays
using StructArrays
using Unitful

using LegendDataManagement

using ProgressMeter: @showprogress
using PropertyFunctions: PropertyFunction, @pf

using Tables: columns

using LegendDataTypes: fast_flatten

include("flatten_over_channels.jl")
include("build_global_events.jl")

end # module
