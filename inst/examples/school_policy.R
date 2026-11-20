library(fasttwfeR)
school_data <- read.csv(system.file("extdata", "school_policy.csv",
                                    package = "fasttwfeR"))
fit <- fasttwfe(school_data, "test_score", "year", "school",
                   gvar = "first_treat", covariates = "baseline_score",
                   data_mode = "rcs", inference = c("native", "jackknife"))
example_result <- data.frame(estimate = fit$estimate, cluster_se = fit$native_se,
                             jackknife_se = fit$jackknife_se, jackknife_p = fit$jackknife_p)
print(example_result, digits = 7, row.names = FALSE)
