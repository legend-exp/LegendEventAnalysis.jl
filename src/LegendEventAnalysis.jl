# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

__precompile__(true)

module LegendEventAnalysis

using ArraysOfArrays
using StructArrays
import Statistics
using Statistics: mean, median, std
using Unitful, UnitfulAtomic
using Dates: DateTime, unix2datetime
import Dates
import Measurements
using Measurements: measurement

using LegendDataManagement

using IntervalSets: AbstractInterval, ClosedInterval, leftendpoint, rightendpoint
using ProgressMeter
using PropertyFunctions: PropertyFunction, @pf, filterby, PropSelFunction

using Tables: columns
using TypedTables: Table
using Unitful: RealOrRealQuantity as RealQuantity
using LegendDataTypes: fast_flatten


const AbstractDataStore = Union{AbstractDict, NamedTuple}

include("flatten_over_channels.jl")
include("build_global_events.jl")
include("calibrate_geds.jl")
include("calibrate_smps.jl")
include("calibrate_pmts.jl")
include("calibrate_aux.jl")
include("calibrate_all.jl")
include("time_evolution.jl")

end # module
