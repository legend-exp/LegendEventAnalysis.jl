## Time evolution

`TimeEvolution(time, values)` holds one quantity sampled at absolute times.
`time` is stored as `DateTime`; unitful unix timestamps (as found in the
`hit` and `evt` tiers) are converted on construction. `values` may carry
units and `Measurement` uncertainties. Timestamps must be sorted.

Three functions derive new series from existing data, and each result is
again a `TimeEvolution`, so they compose:

```julia
using LegendEventAnalysis, Unitful

ts = TimeEvolution(timestamp, e_pulser)        # raw pulser energies, one per event
s = smooth(ts, 1u"hr")                      # mean ± std per hour, one point per window
g = relative(s, 6u"hr")                     # deviation in percent from the first 6 hours
r = rate(timestamp, 10u"minute")            # events per second, √N uncertainties
```

- `smooth(ts, Δt)` averages over consecutive windows of width `Δt` and
  reports the spread of each window as the uncertainty.
- `relative(ts, reference)` normalizes to the mean of the first `reference`
  samples (an `Integer`) or of the samples within a time span after the first
  timestamp. Units cancel; uncertainties propagate.
- `rate(timestamps, Δt)` counts events per window and reports the rate in Hz.
  Only windows fully covered by the data span are included.

Window widths and reference spans accept `Dates.Period`s (`Hour(1)`) as well
as unitful times (`1u"hr"`).

`smooth` takes any statistic as a third argument; `mean` (the default) adds
the window spread as uncertainty and `extrema` yields `ClosedInterval`s, which
plot as an envelope:

```julia
smooth(ts, 1u"hr", median)             # median per window
smooth(ts, 1u"hr", std)                # spread per window, e.g. to monitor the resolution
smooth(ts, 1u"hr", maximum)            # or minimum, or v -> quantile(v, 0.99)
smooth(ts, 1u"hr", extrema)            # min .. max per window
```

Whole-series statistics and selection:

```julia
mean(ts), std(ts), median(ts), extrema(ts), minimum(ts), maximum(ts)
findmax(ts), findmin(ts)               # (value, time) of the extreme sample
ts[t1 .. t2]                           # samples within a ClosedInterval of DateTime or unix times
filter(>(1005u"keV"), ts)              # samples whose values satisfy a predicate
```

## Time-evolution plots

Loading LegendMakie and a Makie backend enables `lplot` and `lplot!` for
`TimeEvolution`. The examples below use

```julia
using LegendMakie, CairoMakie, Dates, Unitful

ts = TimeEvolution(timestamp, e_pulser)    # one pulser energy per event
g = relative(smooth(ts, 1u"hr"), 6u"hr")   # deviation in percent from the first 6 hours
r = rate(timestamp, 30u"minute")           # event rate in Hz
```

### Basic plot

By default the values are drawn with `lines!` and their uncertainties as
bands at `sigmas = (3, 1)`:

```julia
lplot(g)
lplot(g; ylabel = "ΔE (%)", title = "Gain stability", figsize = (1000, 400))
```

### Uncertainty styles

```julia
lplot(g; sigmas = (1,))                                        # a single ±1σ band
lplot(g; sigmas = ())                                          # no bands
lplot(g; uncertainty_style = :bars, plot = scatter!)           # ±1σ error bars
lplot(g; uncertainty_style = :bars, whiskerwidth = 0)          # error bars without caps
lplot(g; uncertainty_style = :none)                            # values only
```

Values without uncertainties (for example, a `TimeEvolution` built from plain
numbers) draw the line or markers only.

### Plot primitive and Makie keywords

`plot` selects the Makie function used for the values; keywords that are not
consumed by `lplot` are forwarded to it:

```julia
lplot(g; plot = scatter!, markersize = 6)
lplot(g; plot = scatterlines!, linestyle = :dash, marker = :diamond)
lplot(g; color = :red, linewidth = 3, label = "Pulser")
```

### Boundaries, axis options, watermark

```julia
lplot(g; boundaries = run_start_times)                         # dashed vertical lines
lplot(g; ylims = (-0.3, 0.3))
lplot(g; axis = (; xticklabelrotation = π/6))                  # any other Axis attribute
lplot(g; watermark = true, watermark_options = (; position = "outer top", final = false))
```

### Overlaying series

`lplot!(ax, ts)` adds a series to an existing axis:

```julia
fig = lplot(g; label = "Pulser", ylabel = "ΔE (%)")
ax = content(fig[1, 1])
lplot!(ax, relative(smooth(cusp, 1u"hr"), 6u"hr"); color = :red, label = "E_cusp")
axislegend(ax)
fig
```

### Envelopes

Interval-valued series from `smooth(ts, Δt, extrema)` draw both edges and the
filled band between them:

```julia
fig = lplot(smooth(ts, 1u"hr", extrema); label = "min–max", ylabel = "Pulser energy")
lplot!(content(fig[1, 1]), smooth(ts, 1u"hr", median); color = :red, label = "median")
axislegend(content(fig[1, 1]))
```

### Several panels

`lplot!(ts; row, col)` adds a new axis to the current figure. Unitful values
label the axis with their unit:

```julia
fig = lplot(g; ylabel = "ΔE (%)")
lplot!(r; row = 2, plot = scatter!, uncertainty_style = :bars, ylabel = "Rate")
lplot!(smooth(ts, 1u"hr"); row = 1, col = 2, ylabel = "Pulser energy")
fig
```
