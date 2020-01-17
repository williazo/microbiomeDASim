set.seed(011520)
id_list <- lapply(seq_len(30), function(i){
    obs <- sample(seq_len(10), size=1)
    id_rep <- rep(i, obs)
})

time_interval <- c(0, 10)
time_list <- lapply(id_list, function(x){
    time_len <- length(x)
    times <- runif(time_len, min=time_interval[1], max=time_interval[2])
    times <- times[order(times)]
})

group_list <- lapply(id_list, function(x){
    group_len <- length(x)
    tx_ind <- sample(seq_len(2), 1)
    tx_group <- ifelse(tx_ind==1, "Control", "Treatment")
    groups <- rep(tx_group, group_len)
})
id <- unlist(id_list)
group <- factor(unlist(group_list), levels = c("Control", "Treatment"))
time <- unlist(time_list)

test_that("missmatch of vector length", {
    id_mm <- id[-1]
    time_mm <- time[-1]
    group_mm <- group[-1]
    expect_error(mvrnorm_sim_obs(id_mm, time, group, ref="Control", control_mean=2,
                                    sigma=1, rho=0.9, corr_str="comp",
                                    func_form="L_up", beta=0.5, IP=5),
                "id, time, and group must all be length N")
    expect_error(mvrnorm_sim_obs(id, time_mm, group, ref="Control", control_mean=2,
                                    sigma=1, rho=0.9, corr_str="comp",
                                    func_form="L_up", beta=0.5, IP=5),
                "id, time, and group must all be length N")
    expect_error(mvrnorm_sim_obs(id, time, group_mm, ref="Control", control_mean=2,
                                    sigma=1, rho=0.9, corr_str="comp",
                                    func_form="L_up", beta=0.5, IP=5),
                "id, time, and group must all be length N")
})

test_that("non-factor group", {
    group_char <- as.character(group)
    expect_error(mvrnorm_sim_obs(id, time, group_char, ref="Control", control_mean=2,
                                    sigma=1, rho=0.9, corr_str="comp",
                                    func_form="L_up", beta=0.5, IP=5),
                "group must be factor variable")
})

test_that("improper ref", {
    expect_error(mvrnorm_sim_obs(id, time, group, ref="Active", control_mean=2,
                                    sigma=1, rho=0.9, corr_str="comp",
                                    func_form="L_up", beta=0.5, IP=5),
                "ref must be specified as one of group values")
})

test_that("non-binary group", {
    tert_val <- sample(seq_len(3), size=length(id), replace=TRUE)
    group_tert <- factor(tert_val, levels=seq_len(3),
                            labels=paste0("group", seq_len(3)))
    expect_error(mvrnorm_sim_obs(id, time, group_tert, ref="group1", control_mean=2,
                                    sigma=1, rho=0.9, corr_str="comp",
                                    func_form="L_up", beta=0.5, IP=5),
                "group must be binary")
})
