# This file is a part of LegendEventAnalysis.jl, licensed under the MIT License (MIT).

using PrecompileTools: @setup_workload, @compile_workload

# Synthetic per-detector data with the concrete types the event-building part of
# `calibrate_all` works on: scalar HPGe columns, SiPM columns holding one
# variable-length vector per trigger, and the `DetectorId` keys of the data store.
@setup_workload begin
    n = 4
    timestamps(offset) = collect((0:n-1) .* 100.0u"μs" .+ offset)
    ged_det = DetectorId.(["V00000A", "V00000B"])
    spm_det = DetectorId.(["S001", "S002"])

    ged_caldata = Dict(det => StructVector(
        timestamp = timestamps(i * 1.0u"μs"),
        e_cusp_ctc_cal = fill(1.0u"keV", n),
        is_physical_trig = trues(n),
    ) for (i, det) in enumerate(ged_det))

    spm_caldata = Dict(det => StructVector(
        timestamp = timestamps(i * 1.0u"μs"),
        trigger_pos = VectorOfVectors([[1.0, 2.0]u"μs" for _ in 1:n]),
        is_dc = VectorOfVectors([[false, true] for _ in 1:n]),
        energy = VectorOfVectors([[1.0, 2.0]u"e_au" for _ in 1:n]),
    ) for (i, det) in enumerate(spm_det))

    @compile_workload begin
        ged_events = build_global_events(ged_caldata, ged_det)
        spm_events = StructArray(map(_fix_vov, columns(build_global_events(spm_caldata, spm_det))))
        events = build_cross_system_events((geds = ged_events, spms = spm_events))
        flag_coincidences(timestamps(0.0u"μs"), timestamps(1.0u"μs"))

        spms = events.spms
        t_wins = fill(ClosedInterval(0.0u"μs", 10.0u"μs"), n)
        _spms_win_sums([:win], _ -> t_wins, "prefix_", spms.trigger_pos, spms.is_dc, spms.energy)
        _lar_cut.(t_wins, spms.trigger_pos, spms.is_dc, spms.energy, Ref(1.0u"e_au"), Ref(1))
        VectorOfVectors(_spm_win_pe_sums.(t_wins, spms.trigger_pos, spms.is_dc, spms.energy))
        _geom_spms_idxs.(Ref(Dict(ged_det[1] => Set(spm_det[1:1]))), fill(ged_det[1:1], n), spms.detector)
    end
end
