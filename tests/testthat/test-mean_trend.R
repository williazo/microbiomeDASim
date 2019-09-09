test_that("functional form misspecification", {
    expect_error(mean_trend(timepoints=seq(0, 10, length.out=20),
                          form="Linear", beta=c(0, -0.5)))

    expect_error(mean_trend(timepoints=seq(0, 10, length.out=20),
                              form="M", beta=c(0, 0.5), IP=c(1)),
                   "IP must have three values")

    expect_error(mean_trend(timepoints=seq(0, 10, length.out=20),
                            form="M", beta=c(0, 0.5), IP=c(0, 1, 2)),
                 "IP cannot equal min or max of timepoints")

    M_test <- mean_trend(timepoints=seq(0, 10, length.out=20),
                         form="M", beta=c(0, 0.5), IP=c(3, 5, 8),
                         plot_trend=FALSE)
    W_test <- mean_trend(timepoints=seq(0, 10, length.out=20),
                          form="W", beta=c(0, -0.1), IP=c(3, 5, 8),
                          plot_trend=FALSE)

    expect_equal(length(M_test), 3)
    expect_equal(names(M_test), c("form", "trend", "beta"))

    expect_error(mean_trend(timepoints=1:5, form="L_down", beta=c(1),
                            IP=3),
                 "For a L_down trend, must specify beta as (beta_yintercept, beta_slope)*")

    expect_warning(mean_trend(timepoints=1:5, form="M", beta=c(0,1),
                            IP=c(3, 2, 4)),
                 "Re-ordering IP from smallest to largest")

    expect_error(mean_trend(timepoints = seq(1, 3,length.out = 100),
                              form = "L_up", beta = 1, IP = c(2, 3, 5),
                              plot_trend = FALSE), "IP must be single value")

    expect_error(mean_trend(timepoints = seq(1, 10,length.out = 100), form = "L_down", beta = c(4, -0.1),
                            IP = NULL, plot_trend = FALSE))

})
