test_that("incomplete mvrnorm arguments", {
    # missing  mvrnorm parameters
    expect_error(gen_norm_microbiome(features=10,
                                diff_abun_features=5))

    comp_ex <- gen_norm_microbiome(features=4,
                              diff_abun_features=2, n_control=5,
                              n_treat=5, control_mean=2,
                              sigma=1, num_timepoints=4, rho=0.8,
                              corr_str="compound", func_form="linear",
                              beta= c(0, 1), missing_pct=0.3,
                              missing_per_subject=2, dis_plot=FALSE)
    #checking that the output is equal to the mrvnorm_object
    expect_equal(length(comp_ex),2)
    #checking the names of the output is equal to mvnorm_object
    expect_equal(names(comp_ex), c("Y", "bug_feat"))
    #checking that differential features is not more than total features specified
    expect_error(names(gen_norm_microbiome(features=1,
                                      diff_abun_features=2, n_control=5,
                                      n_treat=5, control_mean=2, sigma=1,
                                      num_timepoints=4, rho=0.8,
                                      corr_str="compound", func_form="linear",
                                      beta= c(0, 1), missing_pct=0.3,
                                      missing_per_subject=2,
                                      zero_trunc=FALSE)),
                 "differential abundance must be <= to total features")
    exp_neg <-gen_norm_microbiome(features=4,
                             diff_abun_features=2, n_control=10,
                             n_treat=10, control_mean=0, sigma=1,
                             num_timepoints=4, rho=0.8,
                             corr_str="compound", func_form="linear",
                             beta= c(0, 1), missing_pct=0.3,
                             missing_per_subject=2, dis_plot=FALSE,
                             zero_trunc=FALSE)
    #expect all values to have at least one negative value
    expect_true(sum(apply(exp_neg$Y, 1, min) < 0) == 4)

    #all features differentially abundant
    expect_warning(names(gen_norm_microbiome(features=2,
                                        diff_abun_features=2, n_control=20,
                                        n_treat=20, control_mean=2, sigma=1,
                                        num_timepoints=6, rho=0.8,
                                        corr_str="compound", func_form="linear",
                                        beta= c(0, 1), missing_pct=0.3,
                                        missing_per_subject=2,
                                        zero_trunc=FALSE)),
                   "all features will be simulated with differential abundance")

    #no features differentially abundant
    expect_warning(names(gen_norm_microbiome(features=2,
                                        diff_abun_features=0, n_control=20,
                                        n_treat=20, control_mean=2, sigma=1,
                                        num_timepoints=6, rho=0.8,
                                        corr_str="compound", func_form="linear",
                                        beta= c(0, 1), missing_pct=0.3,
                                        missing_per_subject=2,
                                        zero_trunc=FALSE)),
                   "no differential abundance features specified")

    #no features differentially abundant
    expect_error(names(gen_norm_microbiome(features=0,
                                      diff_abun_features=0, n_control=20,
                                      n_treat=20, control_mean=2, sigma=1,
                                      num_timepoints=6, rho=0.8,
                                      corr_str="compound", func_form="linear",
                                      beta= c(0, 1), missing_pct=0.3,
                                      missing_per_subject=2,
                                      zero_trunc=FALSE)),
                 "must specify features greater than zero")


})

