# Two-way fixed-effects DiD in R

Estimate a policy coefficient with school (or state) and year fixed effects,
then compare conventional clustered uncertainty with a delete-cluster jackknife.
The package fits one pooled regression with a treatment indicator and optional
covariates. Use it when that regression is the specification you intend to study.

## Install

Use R 4.1 or later. Estimation uses base R; it does not call Python or Stata.
Install this source version with:

```r
install.packages("remotes")  # once, if needed
remotes::install_github("mattdwebb/fasttwfeR")
```

Alternatively, download the repository, open a terminal in its directory, and run
`R CMD INSTALL .`. For reproducibility, record the installed version
with `packageVersion("fasttwfeR")` and the source commit.

## A school-policy example

The included `school_policy.csv` contains 960 pupil observations: ten pupils
in each of 12 schools over 2010–2017. Four schools adopt in 2013, four in 2015,
and four never adopt. `test_score` is the outcome in score points and
`baseline_score` is a predetermined score covariate. Pupils are repeated cross
sections; their within-school-year numbers are not longitudinal identifiers.

These are fictional teaching data. The generator adds a two-point policy effect
and varying disturbances; the example is not evidence about a real education
policy. [Data definitions and the generator](inst/extdata/DATA_GUIDE.md) are
included, and the same input is used in the Python and Stata examples.

```r
library(fasttwfeR)
school_data <- read.csv(system.file("extdata", "school_policy.csv",
                                    package = "fasttwfeR"))
fit <- fasttwfe(school_data, "test_score", "year", "school",
                   gvar = "first_treat", covariates = "baseline_score",
                   data_mode = "rcs", inference = c("native", "jackknife"))
example_result <- data.frame(estimate = fit$estimate, cluster_se = fit$native_se,
                             jackknife_se = fit$jackknife_se, jackknife_p = fit$jackknife_p)
print(example_result, digits = 7, row.names = FALSE)
```

Actual output, rounded to six decimal places:

| Estimate | Clustered SE | Jackknife SE | Jackknife p-value |
| --- | --- | --- | --- |
| 2.012599 | 0.367339 | 0.401294 | 0.000393 |

The estimated adjusted effect is **2.01 score points**. Its jackknife
standard error is **0.40 points**. The p-value tests a zero effect using a
Student t distribution with 11 degrees of freedom, one less than the 12 schools.
The 960 pupil records do not create 960 independent policy assignments.
The two standard errors describe uncertainty about the same reported estimate.
The clustered calculation uses this estimator's native variance formula. The
jackknife recomputes the estimate after deleting each school, including nuisance
models and aggregation weights.

## Adapt the example

Omit `covariates` to estimate without additional controls. To change the
covariate specification, pass the names of numeric columns, for example
`covariates = c("baseline_score", "class_size")` after adding `class_size` to
your own data. Encode categorical controls before passing them. Covariates should
be chosen to suit the identification assumptions; a variable affected by the
policy is not made appropriate by including it in a regression.

`first_treat` repeats each school's first policy year on every row. Zero means
never treated. Treatment must remain on after adoption. Keep all necessary
pre-policy observations and check that the required control comparisons exist.
Missing or unsupported inputs produce explicit errors instead of silently
changing the specification.

Save the displayed results with `write.csv(example_result, "school_results.csv",
row.names = FALSE)`. The full `fit` object retains the estimate, inference results,
and the estimator-specific diagnostics described in `?fasttwfe`.

## Inference and supported scope

TWFE uses all observations in a pooled regression. With staggered adoption and
heterogeneous effects, its coefficient need not equal a positively weighted
average treatment effect. The jackknife changes the uncertainty calculation;
it does not change the regression's target.

The jackknife uses all input inference clusters. If deleting a school removes
one adoption cohort but leaves an estimable aggregate, the supported effects and
their weights are recomputed. If a required deletion has no estimable aggregate,
the jackknife fails explicitly; failed deletions are never discarded.

## Help, validation and citation

Use `help("fasttwfe", package = "fasttwfeR")` for arguments and return values.
The [method specification](METHOD_SPECIFICATION.md), [validation record](VALIDATION.md),
and [developer reference](docs/DEVELOPER_REFERENCE.md) give additional detail.
These are supported subsets, not promises of every feature in the original
estimators' software. Cite the relevant method and use `citation("fasttwfeR")`
for the software version you used. Original code is distributed under the
[package license](LICENSE); third-party notices retain their own terms.

From a source checkout, run `R CMD INSTALL .` followed by
`Rscript tests/readme.R` to reproduce and check this example. The executable
source is [school_policy.R](inst/examples/school_policy.R).

Current candidate checks: [release validation](RELEASE_VALIDATION.md).
