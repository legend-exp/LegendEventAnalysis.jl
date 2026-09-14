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

`gain_stability(time, energy; ...)` creates a time-series struct for one energy
series. Call it independently for each series that you want to compare. It
computes the normalization, rolling mean and standard deviation before plotting,
preserving the input times.
The caller must remove nonfinite time/energy pairs and provide nonempty data,
a nonzero first energy, and positive `n_ref` and `n_smooth`. `gain_stability`
does not filter or validate these inputs.
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

Uncertainties are optional: `TimeSeriesStruct(time, values)` plots values alone.
For a struct that already has uncertainties, use `uncertainty_style = :none` to
hide them, `:bars` for one-sigma error bars, or `:band` (the default) for bands.
`sigmas` applies to bands; `whiskerwidth` controls error-bar caps.

The A/E mean, width, fast shifts, and slow drifts each fit the same struct.
For example, one panel with run labels can use:

```julia
series = TimeSeriesStruct(run_times, mu; uncertainty = mu_errors)
lplot!(ax, series; plot = scatterlines!, uncertainty_style = :bars,
    color = :blue, marker = :circle, linestyle = :dash)
```

Numeric run positions with custom `xticks` also work. Partition shading,
reference lines, and the four-panel layout remain ordinary Makie operations.
