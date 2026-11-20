\name{fasttwfe}
\alias{fasttwfe}
\title{Static TWFE Estimation and Inference}
\usage{
fasttwfe(data, outcome, time, cluster, treatment = NULL, gvar = NULL,
  covariates = character(), weights = NULL, unit = NULL,
  data_mode = "auto", inference = c("native", "jackknife"),
  assignments = NULL, ri_reps = 0L, seed = 1L, null = 0)
}
\arguments{
\item{data}{Data frame; missing observations fail.}
\item{outcome, time, cluster}{Column names as strings.}
\item{treatment}{Binary absorbing treatment column.}
\item{gvar}{Cluster-constant first-treatment column; zero means never.}
\item{covariates}{Numeric covariate column names.}
\item{weights}{Optional positive observation-weight column.}
\item{unit}{Panel unit column, introducing unit FE where supported.}
\item{data_mode}{auto, panel, or rcs.}
\item{inference}{point, native, jackknife, ri_beta, ri_t_native, , .}
\item{assignments}{Alternative cohort-label matrix, columns sorted by cluster ID.}
\item{ri_reps}{Number of unique alternative assignments.}
\item{seed}{RNG seed.}
\item{null}{Hypothesized scalar; RI currently only supports zero.}
}
\value{List with estimate, target, native and jackknife inference, randomization
p-values, deletions, assignments and diagnostics.}
\description{Independent base R development implementation. See package method
and validation documents for precise supported scope.}
