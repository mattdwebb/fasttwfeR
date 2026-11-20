# Run from the package root, or from tests/ against an installed package.
package <- "fasttwfeR"
env <- if (file.exists("R/estimator.R")) {
  e <- new.env(parent = globalenv())
  sys.source("R/estimator.R", envir = e)
  e
} else asNamespace(package)

# Independent full-string reference: no environment names or hash functions.
reference_draws <- function(a, reps, seed) {
  set.seed(seed)
  seen <- paste(a, collapse = ",")
  rows <- list()
  while (length(rows) < reps) {
    draw <- a[sample.int(length(a))]
    key <- paste(draw, collapse = ",")
    if (!key %in% seen) {
      seen <- c(seen, key)
      rows[[length(rows) + 1L]] <- draw
    }
  }
  do.call(rbind, rows)
}

# The original failure: a serialized assignment longer than 10,000 bytes.
a <- c(3, 3, rep(0, 5998))
long_key <- paste(a, collapse = ",")
stopifnot(nchar(long_key, type = "bytes") > 10000L)
expected <- reference_draws(a, 4L, 19L)
set.seed(291)
before <- .Random.seed
actual <- env$.assignments(a, 4L, 19L)
stopifnot(identical(actual, expected), identical(before, .Random.seed))
rm(".Random.seed", envir = .GlobalEnv)
stopifnot(identical(env$.assignments(a, 4L, 19L), expected),
          !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

# Real collision, duplicates, empty strings, Unicode, and distinct long keys.
seen <- new.env(parent = emptyenv())
stopifnot(identical(env$.ri_key_bucket("abc"), env$.ri_key_bucket("b`d")))
env$.ri_key_add(seen, "abc")
stopifnot(!env$.ri_key_has(seen, "b`d"))
env$.ri_key_add(seen, "b`d")
env$.ri_key_add(seen, "abc")
stopifnot(length(seen[[env$.ri_key_bucket("abc")]]) == 2L)
for (key in c("", "\u00e9\u4e2d", long_key, paste0(long_key, ",3"))) {
  env$.ri_key_add(seen, key)
  stopifnot(env$.ri_key_has(seen, key))
}
stopifnot(!env$.ri_key_has(seen, paste0(long_key, ",missing")),
          all(nchar(ls(seen), type = "bytes") < 100L))

# Force every assignment into one bucket: the sampler must still be exact.
colliding <- new.env(parent = env)
colliding$.ri_key_bucket <- function(key) "one-bucket"
for (name in c(".ri_key_has", ".ri_key_add", ".assignments")) {
  fun <- get(name, envir = env)
  environment(fun) <- colliding
  assign(name, fun, envir = colliding)
}
stopifnot(identical(colliding$.assignments(a, 4L, 19L), expected))

# This design has >100,000 possible assignments and takes the hashed path.
# Compare every estimate, covariance/inference result and randomized statistic
# with the same seeded draws from the independent reference sampler.
d <- expand.grid(s = 1:15, t = 1:4, j = 1:2)
cohorts <- rep(c(0, 3, 4), each = 5L)
d$g <- cohorts[d$s]
d$y <- d$s/10 + d$t/3 + 2*(d$g > 0 & d$t >= d$g) +
  sin(d$s*d$t)/5 + cos(d$j + d$s)/7
draws <- reference_draws(cohorts, 3L, 29L)
methods <- c("native", "jackknife", "ri_beta", "ri_t_native")
set.seed(192)
before <- .Random.seed
generated <- env$.run(d, "y", "t", "s", gvar = "g", inference = methods,
                      ri_reps = 3L, seed = 29L)
stopifnot(identical(before, .Random.seed))
supplied <- env$.run(d, "y", "t", "s", gvar = "g", inference = methods,
                     assignments = draws)
generated$runtime <- supplied$runtime <- NULL
stopifnot(identical(generated, supplied),
          all(is.finite(unlist(generated[c("native_p", "jackknife_p", "ri_beta_p",
                                          "ri_t_native_p")]))))

# Long UTF-8 cluster labels must not alter ordering, seeded draws or inference.
labelled <- d
labelled$s <- paste0(sprintf("s%03d-", d$s), strrep("\u00e9", 6000L))
stopifnot(nchar(labelled$s[1L], type = "bytes") > 10000L)
long_fit <- env$.run(labelled, "y", "t", "s", gvar = "g", inference = methods,
                     ri_reps = 3L, seed = 29L)
long_fit$runtime <- NULL
long_fit$cluster_ids <- generated$cluster_ids <- NULL
stopifnot(identical(long_fit, generated))
cat("PASS:", package, "many-cluster keys, collisions, long labels, RNG, seeded estimates and inference\n")
