test_that("linear_trend beta misspecification", {
  expect_error(mvrnorm_sim(n_control=10, n_treat=10, control_mean=2, sigma=1,
                           rho=0.6, num_timepoints=10, t_interval=c(0, 9),
                           func_form="linear", corr_str="ar1",
                           beta=c(0, -0.5, 1), missing_pct=0,
                           missing_per_subject=0),
               "For a linear trend, must specify beta as \\(beta_0, beta_1\\).")

  expect_error(mvrnorm_sim(n_control=10, n_treat=10, control_mean=2, sigma=1,
                           rho=0.6, num_timepoints=10, t_interval=c(0, 9),
                           func_form="linear", corr_str="ar1", beta=c(0),
                           missing_pct=0, missing_per_subject=0),
               "For a linear trend, must specify beta as \\(beta_0, beta_1\\).")
})

test_that("rho misspecification", {
  expect_error(mvrnorm_sim(n_control=10, n_treat=10, control_mean=2,
                           sigma=1, rho="1", num_timepoints=10,
                           t_interval=c(0, 9), func_form="linear",
                           corr_str="ar1", beta=c(0, 1), missing_pct=0,
                           missing_per_subject=0),
               "rho must be numeric")

  expect_error(mvrnorm_sim(n_control=10, n_treat=10, control_mean=2,
                           sigma=1,rho=1.4, num_timepoints=10,
                           t_interval=c(0, 9), func_form="linear",
                           corr_str="ar1", beta=c(0, 1),
                           missing_pct=0, missing_per_subject=0),
               "rho must be in \\[0, 1\\]")
})

test_that("zero truncation", {
    Y_trunc <- mvrnorm_sim(n_control=10, n_treat=10, control_mean=3,
                sigma=2, rho=0.8, num_timepoints=10, t_interval=c(0, 9),
                func_form="linear", corr_str="ar1", beta=c(0, 1),
                missing_pct=0, missing_per_subject=0, zero_trunc=TRUE)
    Y_no_trunc <- mvrnorm_sim(n_control=10, n_treat=10, control_mean=0,
                  sigma=2, rho=0.8, num_timepoints=10, t_interval=c(0, 9),
                  func_form="linear",corr_str="ar1", beta=c(0, 1),
                  missing_pct=0, missing_per_subject=0, zero_trunc=FALSE)
    expect_true(min(Y_trunc$Y)>=0)
    expect_true(min(Y_no_trunc$Y)<0)
})

test_that("number timepoints misspecification", {
  expect_error(mvrnorm_sim(n_control=10, n_treat=10, control_mean=3,
                        sigma=2, rho=0.8, num_timepoints="3", t_interval=1,
                        func_form="linear", corr_str="ar1", beta=c(0, 1),
                        missing_pct=0, missing_per_subject=0, zero_trunc=TRUE),
               "num_timepoints and/or t_interval not numeric")
  expect_error(mvrnorm_sim(n_control=10, n_treat=10, control_mean=3,
                        sigma=2, rho=0.8, num_timepoints=10,
                        t_interval="5",func_form="linear", corr_str="ar1",
                        beta=c(0, 1), missing_pct=0, missing_per_subject=0,
                        zero_trunc=TRUE),
               "num_timepoints and/or t_interval not numeric")
})


test_that("time interval misspecification", {
  expect_error(mvrnorm_sim(n_control=10, n_treat=10, control_mean=3,
                        sigma=2, rho=0.8, num_timepoints=10, t_interval=1,
                        func_form="linear", corr_str="ar1", beta=c(0, 1),
                        missing_pct=0, missing_per_subject=0, zero_trunc=TRUE),
  "time interval must be numeric vector of length 2")
  expect_error(mvrnorm_sim(n_control=10, n_treat=10, control_mean=3,
                        sigma=2, rho=0.8, num_timepoints=10,
                        t_interval=c(0, Inf),func_form="linear", corr_str="ar1",
                        beta=c(0, 1), missing_pct=0, missing_per_subject=0,
                        zero_trunc=TRUE),
               "time interval must be finite")
})

