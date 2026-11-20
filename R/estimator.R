# Exact set membership for RI assignment keys. R environment names are
# symbols limited to 10000 bytes; the full key must be a value, not a name.
# Short hash buckets keep lookups fast. Comparing full keys within each
# bucket prevents hash collisions from discarding distinct assignments.
.ri_key_bucket <- function(key) {
  code <- as.numeric(utf8ToInt(enc2utf8(key)))
  n <- length(code)
  index <- seq_len(n)
  forward <- sum((code * (index %% 65521 + 1)) %% 2147483647) %% 2147483647
  reverse <- sum((code * (rev(index) %% 65519 + 7)) %% 2147483629) %% 2147483629
  sprintf("h%.0f:%.0f:%d", forward, reverse, n)
}

.ri_key_has <- function(seen, key) {
  key %in% seen[[.ri_key_bucket(key)]]
}

.ri_key_add <- function(seen, key) {
  bucket <- .ri_key_bucket(key)
  keys <- seen[[bucket]]
  if (!key %in% keys) seen[[bucket]] <- c(keys, key)
  invisible(NULL)
}

.family <- "twfe"

.fail <- function(reason) stop(reason, call. = FALSE)

.absorb <- function(z, groups, w, tol = 1e-12, maxit = 10000L) {
  z <- as.matrix(z)
  scale <- pmax(1, apply(abs(z), 2L, max))
  for (iteration in seq_len(maxit)) {
    old <- z
    for (g in groups) {
      mass <- as.vector(rowsum(w, g, reorder = FALSE))
      lev <- unique(g)
      means <- rowsum(z * w, g, reorder = FALSE) / mass
      z <- z - means[match(g, lev), , drop = FALSE]
    }
    if (max(sweep(abs(z - old), 2L, scale, "/")) < tol)
      return(list(residual = z, iterations = iteration))
  }
  .fail("NUMERICAL_FAILURE: FE absorption did not converge")
}

.solve_fit <- function(x, y, w, allow_rank = FALSE) {
  xw <- x * sqrt(w)
  s <- svd(xw)
  active <- s$d > max(dim(xw)) * .Machine$double.eps * max(s$d)
  if (!allow_rank && sum(active) != ncol(x))
    .fail("IDENTIFICATION_FAILURE: rank deficient design")
  if (!any(active)) .fail("IDENTIFICATION_FAILURE: zero rank")
  pinv <- (s$v[, active, drop = FALSE] / rep(s$d[active], each = ncol(x))) %*%
    t(s$u[, active, drop = FALSE])
  list(b = as.vector(pinv %*% (sqrt(w) * y)), pinv = pinv, rank = sum(active))
}

.dummies <- function(g, drop = FALSE) {
  k <- max(g)
  z <- matrix(0, length(g), k)
  z[cbind(seq_along(g), g)] <- 1
  if (drop) z <- z[, -1L, drop = FALSE]
  z
}

.prepare <- function(data, outcome, time, cluster, treatment, gvar, covariates,
                     weights, unit, data_mode) {
  columns <- unique(c(outcome, time, cluster, treatment, gvar, covariates, weights, unit))
  if (!is.data.frame(data) || !all(columns %in% names(data)))
    .fail("INVALID_INPUT: data frame or column names")
  if (nrow(data) == 0L || anyNA(data[, columns, drop = FALSE]))
    .fail("INVALID_INPUT: missing data must be handled explicitly")
  numeric_columns <- unique(c(outcome, time, treatment, gvar, covariates, weights))
  if (!all(vapply(data[numeric_columns], is.numeric, logical(1))) ||
      !all(is.finite(as.matrix(data[numeric_columns]))))
    .fail("INVALID_INPUT: nonfinite or nonnumeric model variables")
  if (!data_mode %in% c("auto", "panel", "rcs")) .fail("UNSUPPORTED: data_mode")
  if (data_mode == "auto") data_mode <- if (is.null(unit)) "rcs" else "panel"
  if (data_mode == "panel" && is.null(unit)) .fail("INVALID_INPUT: panel requires unit")
  if (data_mode == "panel" && anyDuplicated(data[c(unit, time)]))
    .fail("INVALID_INPUT: duplicate unit-time")
  ids <- sort(unique(data[[cluster]]))
  c <- match(data[[cluster]], ids)
  t <- data[[time]]
  if (is.null(treatment) && is.null(gvar)) .fail("INVALID_INPUT: supply treatment or gvar")
  if (!is.null(treatment) && !all(data[[treatment]] %in% c(0, 1)))
    .fail("INVALID_INPUT: treatment must be binary")
  a <- numeric(length(ids))
  for (i in seq_along(ids)) {
    ix <- c == i
    if (!is.null(gvar)) {
      values <- unique(data[[gvar]][ix])
      if (length(values) != 1L || values < 0) .fail("INVALID_INPUT: cluster-constant gvar required")
      a[i] <- values
    } else {
      a[i] <- if (any(data[[treatment]][ix] == 1)) min(t[ix & data[[treatment]] == 1]) else 0
    }
    d <- as.numeric(a[i] > 0 & t[ix] >= a[i])
    if (!is.null(treatment) && any(d != data[[treatment]][ix]))
      .fail("INVALID_INPUT: inconsistent or nonabsorbing treatment")
    if (a[i] > 0 && (!any(t[ix] < a[i]) || !any(t[ix] >= a[i])))
      .fail("IDENTIFICATION_FAILURE: each treated cluster requires pre and post")
  }
  if (!any(a > 0)) .fail("IDENTIFICATION_FAILURE: no treatment")
  w <- if (is.null(weights)) rep(1, nrow(data)) else data[[weights]]
  if (any(w <= 0)) .fail("INVALID_INPUT: weights must be positive")
  fe <- if (is.null(unit) || data_mode == "rcs") c else match(data[[unit]], sort(unique(data[[unit]])))
  if (data_mode == "panel" && any(vapply(split(c, fe), function(z) length(unique(z)) != 1L, logical(1))))
    .fail("INVALID_INPUT: units must nest in clusters")
  list(y = data[[outcome]], t = t, c = c, fe = fe, a = a,
       x = as.matrix(data[covariates]), w = w, ids = ids, mode = data_mode)
}

.fit <- function(p, a, omit = NULL, native = FALSE) {
  keep <- if (is.null(omit)) rep(TRUE, length(p$y)) else p$c != omit
  y <- p$y[keep]; w <- p$w[keep]; x <- p$x[keep, , drop = FALSE]
  c <- match(p$c[keep], unique(p$c[keep]))
  f <- match(p$fe[keep], unique(p$fe[keep]))
  t <- match(p$t[keep], sort(unique(p$t[keep])))
  g <- a[p$c[keep]]
  d <- as.numeric(g > 0 & p$t[keep] >= g)
  n <- length(y); G <- length(unique(c))
  if (!any(d == 1) || !any(d == 0)) .fail("IDENTIFICATION_FAILURE: treatment support")
  sums <- function(z) as.vector(rowsum(z, c, reorder = FALSE))
  multiplier <- 1
  scores <- NULL
  iterations <- 0L
  if (.family == "twfe") {
    z <- .absorb(cbind(y, x, d), list(f, t), w)
    iterations <- z$iterations
    yy <- z$residual[, 1L]
    xx <- z$residual[, -1L, drop = FALSE]
    fit <- .solve_fit(xx, yy, w)
    theta <- tail(fit$b, 1)
    if (native) {
      # Rank of the additive FE design via connected components, without dummies.
      parent <- seq_len(max(f) + max(t))
      root <- function(i) { while (parent[i] != i) i <- parent[i]; i }
      for (i in seq_len(n)) {
        u <- root(f[i]); v <- root(max(f) + t[i])
        if (u != v) parent[v] <- u
      }
      components <- length(unique(vapply(seq_along(parent), root, integer(1))))
      rank <- max(f) + max(t) - components + ncol(xx)
      if (G < 2 || n <= rank) .fail("IDENTIFICATION_FAILURE: residual degrees of freedom")
      scores <- sums(sqrt(w) * (yy - as.vector(xx %*% fit$b)) * as.vector(tail(fit$pinv, 1L)))
      multiplier <- G/(G-1) * (n-1)/(n-rank)
    }
  } else {
    q <- ncol(x) + 1L
    width <- if (.family == "flex") (max(f) + max(t) - 1L + sum(d)) * q else max(f) + max(t) - 1L + ncol(x)
    if (n * width > 12000000) .fail("UNSUPPORTED: bounded reference backend exceeds 12 million entries")
    base <- cbind(.dummies(f), .dummies(t, TRUE))
    if (.family %in% c("did2s", "bjs")) {
      z <- cbind(base, x)
      u <- d == 0
      fit <- .solve_fit(z[u, , drop = FALSE], y[u], w[u])
      e <- y - as.vector(z %*% fit$b)
      mass <- sum(w * d)
      theta <- sum(w * d * e) / mass
      if (native) {
        h <- colSums(z * (w * d)) / mass
        direction <- as.vector(t(fit$pinv) %*% h) * sqrt(w[u])
        if (.family == "did2s") {
          score <- w * d * (e - theta) / mass
          score[u] <- -direction * e[u]
        } else {
          if (any(w != w[1])) .fail("UNSUPPORTED: BJS native nonuniform weights")
          residual <- e
          for (gg in unique(g[d == 1])) for (tt in unique(t[g == gg & d == 1])) {
            ix <- g == gg & t == tt & d == 1
            counts <- tabulate(c[ix], nbins = G)
            smart <- counts[c[ix]]
            residual[ix] <- e[ix] - sum(smart * e[ix]) / sum(smart)
          }
          score <- w * d * residual / mass
          score[u] <- -direction * residual[u]
        }
        scores <- sums(score)
      }
    } else {
      if (max(f) != G) .fail("UNSUPPORTED: FLEX requires cluster fixed effects (omit unit or use rcs)")
      v <- cbind(1, x)
      pieces <- lapply(seq_len(ncol(base)), function(j) base[, j] * v)
      cell <- interaction(f, t, drop = TRUE)
      cells <- unique(cell[d == 1])
      post <- lapply(cells, function(cc) as.numeric(cell == cc) * v)
      z <- do.call(cbind, c(pieces, post))
      fit <- .solve_fit(z, y, w, TRUE)
      offset <- ncol(base) * q
      contrast <- numeric(ncol(z))
      effect <- numeric(n)
      mass <- sum(w * d)
      for (j in seq_along(cells)) {
        ix <- cell == cells[j]
        at <- offset + (j - 1L) * q + seq_len(q)
        contrast[at] <- sum(w[ix]) / mass * colMeans(v[ix, , drop = FALSE])
        effect[ix] <- as.vector(v[ix, , drop = FALSE] %*% fit$b[at])
      }
      # Estimability is required even when post-cell covariate columns are redundant.
      if (max(abs(as.vector(contrast %*% fit$pinv %*% (z * sqrt(w))) - contrast)) > 1e-8)
        .fail("IDENTIFICATION_FAILURE: FLEX contrast is not estimable")
      theta <- sum(contrast * fit$b)
      if (native) {
        if (any(w != w[1])) .fail("UNSUPPORTED: FLEX native nonuniform weights")
        if (G < 2 || n <= fit$rank) .fail("IDENTIFICATION_FAILURE: residual degrees of freedom")
        influence <- as.vector(contrast %*% fit$pinv)
        scores <- sums(sqrt(w) * (y - as.vector(z %*% fit$b)) * influence +
                         w * d * (effect - theta) / mass)
        multiplier <- G/(G-1) * (n-1)/(n-fit$rank)
      }
    }
  }
  se <- if (native) sqrt(sum(scores^2) * multiplier) else NA_real_
  if (!is.finite(theta)) .fail("NUMERICAL_FAILURE: nonfinite estimate")
  list(estimate = theta, se = se, scores = scores, multiplier = multiplier,
       iterations = iterations)
}

.jack <- function(p, a) {
  G <- length(a)
  if (G < 2) .fail("IDENTIFICATION_FAILURE: jackknife needs at least two clusters")
  deletions <- vapply(seq_len(G), function(i) .fit(p, a, i)$estimate, numeric(1))
  se <- sqrt((G-1)/G * sum((deletions - mean(deletions))^2))
  if (!is.finite(se) || se <= 64 * .Machine$double.eps * max(1, max(abs(deletions))))
    .fail("NUMERICAL_FAILURE: zero jackknife variance")
  list(se = se, deletions = deletions)
}

.assignments <- function(a, reps, seed) {
  if (length(reps) != 1L || !is.finite(reps) || reps < 1 || reps != floor(reps))
    .fail("INVALID_INPUT: ri_reps must be a positive integer")
  counts <- table(a)
  logspace <- lgamma(length(a)+1) - sum(lgamma(as.numeric(counts)+1))
  space <- exp(logspace)
  key <- function(x) paste(x, collapse = ",")
  observed <- key(a)
  if (space <= 100000) {
    values <- as.numeric(names(counts))
    rows <- list()
    visit <- function(prefix, left) {
      if (sum(left) == 0) {
        if (key(prefix) != observed) rows[[length(rows)+1L]] <<- prefix
        return(invisible(NULL))
      }
      for (j in which(left > 0)) {
        nextleft <- left; nextleft[j] <- nextleft[j] - 1L
        visit(c(prefix, values[j]), nextleft)
      }
    }
    visit(numeric(), as.numeric(counts))
    all <- do.call(rbind, rows)
    if (is.null(all)) .fail("IDENTIFICATION_FAILURE: constant assignment")
    if (reps >= nrow(all)) return(all)
  } else all <- NULL
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  if (!is.null(all)) return(all[sample.int(nrow(all), reps), , drop = FALSE])
  if (reps > 100000) .fail("UNSUPPORTED: at most 100000 random assignments")
  seen <- new.env(hash = TRUE, parent = emptyenv())
  .ri_key_add(seen, observed)
  rows <- vector("list", reps); j <- 0L; tries <- 0L
  while (j < reps) {
    tries <- tries + 1L
    if (tries > 100 * reps) .fail("NUMERICAL_FAILURE: assignment sampler exhausted")
    draw <- a[sample.int(length(a))]
    k <- key(draw)
    if (!.ri_key_has(seen, k)) {
      .ri_key_add(seen, k); j <- j + 1L; rows[[j]] <- draw
    }
  }
  do.call(rbind, rows)
}

.run <- function(data, outcome, time, cluster, treatment = NULL, gvar = NULL,
                 covariates = character(), weights = NULL, unit = NULL,
                 data_mode = "auto", inference = c("native", "jackknife"),
                 assignments = NULL, ri_reps = 0L, seed = 1L, null = 0) {
  started <- proc.time()[3L]
  allowed <- c("point", "native", "jackknife", "ri_beta", "ri_t_native")
  if (!all(inference %in% allowed)) .fail("UNSUPPORTED: inference method")
  if (length(null) != 1L || !is.finite(null)) .fail("INVALID_INPUT: null")
  if (null != 0 && any(grepl("^ri", inference)))
    .fail("UNSUPPORTED: randomization currently tests the zero sharp null only")
  p <- .prepare(data, outcome, time, cluster, treatment, gvar, covariates, weights, unit, data_mode)
  G <- length(p$a)
  native <- any(c("native", "ri_t_native") %in% inference)
  jk <- any(c("jackknife") %in% inference)
  fit <- .fit(p, p$a, native = native)
  result <- list(estimate = fit$estimate, native_se = NA_real_, native_stat = NA_real_,
                 native_p = NA_real_, native_df = NA_real_,
                 jackknife_se = NA_real_, jackknife_t = NA_real_, jackknife_p = NA_real_,
                 jackknife_df = NA_real_, ri_beta_p = NA_real_, ri_t_native_p = NA_real_,
                 deletions = NULL, assignments = NULL,
                 nobs = length(p$y), nclusters = G, cluster_ids = p$ids,
                 data_mode = p$mode, covariates_used = covariates,
                 target_label = if (.family == "twfe") "TWFE coefficient" else
                   if (.family == "flex") "Historical state-group FLEX ATT" else "Treated-observation weighted ATT",
                 status = "PASS", estimator = .family, absorption_iterations = fit$iterations)
  if (native) {
    if (fit$se <= 0 || !is.finite(fit$se)) .fail("NUMERICAL_FAILURE: zero native variance")
    result$native_se <- fit$se
    result$native_stat <- (fit$estimate-null)/fit$se
    result$native_df <- if (.family %in% c("twfe", "flex")) G-1 else Inf
    result$native_p <- 2 * pt(-abs(result$native_stat), result$native_df)
    result$cluster_scores <- fit$scores
    result$finite_sample_multiplier <- fit$multiplier
  }
  if (jk) {
    j <- .jack(p, p$a)
    result$jackknife_se <- j$se
    result$deletions <- j$deletions
    result$jackknife_df <- G-1
    result$jackknife_t <- (fit$estimate-null)/j$se
    result$jackknife_p <- 2 * pt(-abs(result$jackknife_t), G-1)
  }
  if (any(c("ri_beta", "ri_t_native") %in% inference)) {
    if (is.null(assignments)) assignments <- .assignments(p$a, ri_reps, seed)
    assignments <- as.matrix(assignments)
    if (!is.numeric(assignments) || ncol(assignments) != G || nrow(assignments) == 0 ||
        any(!is.finite(assignments))) .fail("INVALID_INPUT: assignment matrix")
    if (anyDuplicated(as.data.frame(assignments))) .fail("INVALID_INPUT: duplicate assignment")
    for (i in seq_len(nrow(assignments))) {
      a <- assignments[i, ]
      if (all(a == p$a)) .fail("INVALID_INPUT: observed assignment must be excluded")
      if (!all(sort(a) == sort(p$a))) .fail("INVALID_INPUT: cohort counts must be preserved")
    }
    values <- matrix(NA_real_, nrow(assignments), 2L,
                     dimnames = list(NULL, c("ri_beta", "ri_t_native")))
    for (i in seq_len(nrow(assignments))) {
      a <- assignments[i, ]
      f <- .fit(p, a, native = "ri_t_native" %in% inference)
      values[i, 1L] <- f$estimate
      if ("ri_t_native" %in% inference) {
        if (!is.finite(f$se) || f$se <= 0) .fail("NUMERICAL_FAILURE: randomized native variance")
        values[i, 2L] <- f$estimate/f$se
      }
    }
    obs <- c(fit$estimate, result$native_stat)
    for (j in seq_len(2L)) if (colnames(values)[j] %in% inference) {
      v <- abs(values[, j]); o <- abs(obs[j])
      hits <- sum(v >= o - 64 * .Machine$double.eps * pmax(1, v, o))
      result[[paste0(colnames(values)[j], "_p")]] <- (1+hits)/(nrow(values)+1)
    }
    result$assignments <- assignments
    result$randomized_statistics <- values
  }
  result$runtime <- unname(proc.time()[3L] - started)
  result
}

fasttwfe <- .run
