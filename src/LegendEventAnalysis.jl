# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

__precompile__(true)

module LegendEventAnalysis

using ArraysOfArrays
using StructArrays
using Unitful

using LegendDataManagement

using IntervalSets: AbstractInterval, ClosedInterval, leftendpoint, rightendpoint
using ProgressMeter: @showprogress
using PropertyFunctions: PropertyFunction, @pf, filterby

using Tables: columns

using LegendDataTypes: fast_flatten

include("flatten_over_channels.jl")
include("build_global_events.jl")
include("calibrate_geds.jl")
include("calibrate_smps.jl")
include("calibrate_puls.jl")
include("calibrate_all.jl")

end # module
