# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using LegendEventAnalysis

using LegendTestData
using SolidStateDetectors

# Doctest setup
DocMeta.setdocmeta!(
    LegendEventAnalysis,
    :DocTestSetup,
    :(using LegendEventAnalysis);
    recursive=true,
)

makedocs(
    sitename = "LegendEventAnalysis",
    modules = [LegendEventAnalysis],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://legend-exp.github.io/LegendEventAnalysis.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Extensions" => "extensions.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/legend-exp/LegendEventAnalysis.jl.git",
    forcepush = true,
    push_preview = true,
)
