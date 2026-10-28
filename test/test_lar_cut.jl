# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

using LegendEventAnalysis
using Test

using IntervalSets: ClosedInterval
using StructArrays: StructArray
using Unitful

import LegendEventAnalysis as LEA
using LegendEventAnalysis: DetectorId, Fill

# photoelectrons, `"e"` in the config
const pe = LegendEventAnalysis.LegendDataManagement.UnitfulAtomic.e_au

@testset "LAr cut" begin
    # two HPGe detectors in strings 1 and 2, three SiPMs in the outer and one in the inner barrel
    geds = StructArray(detector = DetectorId.([:V00001A, :V00002A]), detstring = [1, 2], angle = [0.0, 90.0]u"°")
    spms = StructArray(detector = DetectorId.([:S001, :S002, :S003, :S004]), barrel = [:outer, :outer, :outer, :inner],
                       angle = [0.0, 90.0, 180.0, 0.0]u"°")
    dphi_all = (inner = 180.0u"°", outer = 180.0u"°")
    dphi_geom = (inner = 180.0u"°", outer = 45.0u"°")

    @testset "geometric selection" begin
        @test !LEA._is_geometric((dphi_max = dphi_all,))
        @test LEA._is_geometric((dphi_max = dphi_geom,))

        geom = LEA._geom_spms(geds, spms, dphi_geom)
        @test geom[DetectorId(:V00001A)] == Set(DetectorId.([:S001, :S004]))
        @test geom[DetectorId(:V00002A)] == Set(DetectorId.([:S002, :S004]))
        @test LEA._geom_spms(geds, spms, dphi_all)[DetectorId(:V00001A)] == Set(spms.detector)

        bad_spms = StructArray(detector = spms.detector, barrel = spms.barrel, angle = [0.0, NaN, 180.0, 0.0]u"°")
        @test_throws ArgumentError LEA._geom_spms(geds, bad_spms, dphi_geom)
    end

    # two events with one trigger in every SiPM: 1 PE at 1 µs, not discharged
    spm_det = [spms.detector, spms.detector]
    spm_t = [[[1.0u"µs"] for _ in 1:4] for _ in 1:2]
    spm_dc = [[[false] for _ in 1:4] for _ in 1:2]
    spm_pe = [[[1.0pe] for _ in 1:4] for _ in 1:2]
    t_win = ClosedInterval(0.0u"µs", 5.0u"µs")
    geom = LEA._geom_spms(geds, spms, dphi_geom)

    @testset "triggered HPGe detector" begin
        idxs = LEA._geom_spms_idxs.(Ref(geom), [[DetectorId(:V00001A)], [DetectorId(:V00002A)]], spm_det)
        @test idxs == [[1, 4], [2, 4]]
        sums = LEA._spms_win_sum.(Ref((:pe_sum, :multiplicity)), Ref(t_win), view.(spm_t, idxs), view.(spm_dc, idxs), view.(spm_pe, idxs))
        @test getproperty.(sums, :pe_sum) == [2.0pe, 2.0pe]
        @test getproperty.(sums, :multiplicity) == [2, 2]
        @test all(LEA._lar_cut.(Ref(t_win), view.(spm_t, idxs), view.(spm_dc, idxs), view.(spm_pe, idxs), 1.0pe, 1))
    end

    @testset "no triggered HPGe detector" begin
        # forced triggers: the geometric cut sees no SiPMs, so zero sums and no cut, with the same column types
        no_trig = Fill(DetectorId[], 2)
        idxs = LEA._geom_spms_idxs.(Ref(geom), no_trig, spm_det)
        @test all(isempty, idxs)
        sums = LEA._spms_win_sum.(Ref((:pe_sum, :multiplicity)), Ref(t_win), view.(spm_t, idxs), view.(spm_dc, idxs), view.(spm_pe, idxs))
        @test eltype(sums) == NamedTuple{(:pe_sum, :multiplicity), Tuple{typeof(1.0pe), Int}}
        @test getproperty.(sums, :pe_sum) == [0.0pe, 0.0pe]
        @test getproperty.(sums, :multiplicity) == [0, 0]
        @test !any(LEA._lar_cut.(Ref(t_win), view.(spm_t, idxs), view.(spm_dc, idxs), view.(spm_pe, idxs), 1.0pe, 1))

        # one detector of a string stands for the string
        string_idxs = LEA._geom_spms_idxs.(Ref(geom), Fill([DetectorId(:V00002A)], 2), spm_det)
        @test string_idxs == [[2, 4], [2, 4]]
    end
end
