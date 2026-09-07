## Time-series plots

Loading LegendMakie and a Makie backend enables `lplot` for `TimeSeriesStruct`.
Each report contains one quantity versus time, with optional one-sigma
uncertainties. Keep absolute timestamps (for example, `DateTime` values); labels
and plot styling are keyword arguments, separate from the data.

```julia
using LegendEventAnalysis, LegendMakie, CairoMakie

pulser = gain_stability(time, e_pulser; n_ref = 500, n_smooth = 201)
cusp = gain_stability(time, e_cusp; n_ref = 500, n_smooth = 201)
fig = lplot(pulser; label = "Pulser", ylabel = "ΔE (keV)",
    title = "Gain stability", boundaries = run_start_times, figsize = (1000, 450))
ax = content(fig[1, 1])
lplot!(ax, cusp; color = :red, label = "E_cusp")
axislegend(ax)
```

This replaces `lgainstability(time, e_cusp, e_pulser; ...)`
with two independent reports. `gain_stability` computes the normalization,
rolling mean and standard deviation before plotting, preserving the input times.
Nonfinite time/energy pairs are removed separately for each curve; an empty
selection or zero reference energy is rejected.
`Qbb` defaults to `2039.0` (keV); it can also carry units, e.g. `2039u"keV"`.

Use `TimeSeriesStruct(time, values)` directly for absolute values or prepared
averages, and `plot = scatter!` for sparse data such as one average per run:

```julia
report = TimeSeriesStruct(average_time, averages; uncertainty = standard_deviations)
fig = lplot(report; plot = scatter!, ylabel = "Baseline (ADC)", sigmas = ())
```

The default is `plot = lines!`. Uncertainties draw bands at `sigmas = (3, 1)`;
use `(1,)` for only one sigma or `()` to hide them. `boundaries` supplies times
for vertical run/period markers. Axis options go in `axis = (; limits = ...)`;
remaining keywords go to the line or scatter plot. Use `lplot!(ax, report)` to
overlay another quantity, or `lplot!(report; row = 1, col = 2)` for another panel.
Set `watermark = true` to add LEGEND watermarks, with options such as
`watermark_options = (; position = "outer top", final = false)`.
