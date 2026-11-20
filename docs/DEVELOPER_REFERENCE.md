# Earlier development notes

These are the previous landing-page notes, preserved on 19 November 2026 for
provenance. They may describe superseded restrictions. Use the [current README](../README.md),
[inference guide](../INFERENCE.md), and [change history](../CHANGELOG.md) for current behavior.

# fasttwfeR

Native base-R static TWFE; version 0.1.0. Initial implementation milestone,
not full feature parity or official-software certification.

Install: `R CMD INSTALL .`. Runtime: R only, no Python/Stata.

```r
library(fasttwfeR)
d <- read.csv("tests/fixtures/data.csv")
r <- fasttwfe(d, "y", "time", "state", gvar="g", covariates="x",
  inference=c("native","jackknife"))
r[c("estimate","native_se","native_p","jackknife_se","jackknife_p")]
```

Column names are strings. Supply treatment or gvar (0=never), or both for validation.
Treatment must be absorbing and cluster-time constant. Missing/nonfinite data fail.
Panel requires unit with unique unit-time and uses unit FE. Auto without unit
selects RCS with cluster FE. FLEX requires historical cluster-FE grouping.
Positive weights enter WLS and target masses. Covariates enter contemporaneously,
not frozen: additive common slopes except FLEX's interacted slopes. Absorbed
constant covariates explicitly fail rank checks, not silently dropped. Causal
validity of post-treatment covariates is not guaranteed by regression adjustment.

Result list: estimate, target_label, native_* and jackknife_* fields, all deletions,
RI p-values and randomized statistics, assignments, counts, cluster IDs and runtime.
Unrequested inference is NA; failures raise classified errors for harness isolation.

Unsupported: event studies/windows, distinct assignment units, harmonized targets,
Subtract/Add placebo procedures, anticipation, alternate cohort grouping/aggregation.
Non-TWFE backends are bounded dense references, not large-N optimized ports.
Test: `Rscript tests/run.R .`. Frozen Python is development-only.
