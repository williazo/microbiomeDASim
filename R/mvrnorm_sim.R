#' Simulate Microbiome Longitudinal Data from Multivariate Random Normal
#'
#' This function is used in the
#' \code{\link[microbiomeDASim]{gen_norm_microbiome}} call when the user
#' specified the method as mvrnorm.
#'
#' @param n_control integer value specifying the number of control individuals
#' @param n_treat integer value specifying the number of treated individuals
#' @param control_mean numeric value specifying the mean value for control
#' subjects. all control subjects are assummed to have the same population mean
#'  value.
#' @param sigma numeric value specifying the global population standard
#'  deviation for both control and treated individuals.
#' @param num_timepoints either an integer value specifying the number of
#'  timepoints per subject or a vector of timepoints for each subject. If
#'  supplying a vector the lenght of the vector must equal the total number of
#'  subjects.
#' @param t_interval numeric vector of length two specifying the interval of
#' time from which to draw observatoins \[t_1, t_q\]. Assumed to be equally
#' spaced over the interval unless \code{asynch_time} is set to TRUE.
#' @param rho value for the correlation parameter. must be between \[0, 1\].
#'  see \code{\link[microbiomeDASim]{mvrnorm_corr_gen}} for details.
#' @param corr_str correlation structure selected. see
#'  \code{\link[microbiomeDASim]{mvrnorm_corr_gen}} for details.
#' @param func_form character value specifying the functional form for the
#'  longitduinal mean trend. see \code{\link[microbiomeDASim]{mean_trend}}
#'   for details.
#' @param beta vector value specifying the parameters for the differential
#'     abundance function. see \code{\link[microbiomeDASim]{mean_trend}} for
#'     details.
#' @param IP vector specifying any inflection points. depends on the type of
#' functional form specified. see \code{\link[microbiomeDASim]{mean_trend}} for
#'  details. by default this is set to NULL.
#' @param missing_pct numeric value that must be between \[0, \1] that specifies
#'     what percentage of the individuals will have missing values.
#' @param missing_per_subject integer value specifying how many observations per
#'     subject should be dropped. note that we assume that all individuals must
#'     have baseline value, meaning that the maximum number of
#'  \code{missing_per_subject} is equal to \code{num_timepoints} - 1.
#' @param miss_val value used to induce missingness from the simulated data.
#' by default missing values are assummed to be NA but other common choices
#' include 0.
#' @param dis_plot logical argument on whether to plot the simulated data or
#'     not. by default plotting is turned off.
#' @param plot_trend specifies whether to plot the true mean trend. see
#' \code{\link[microbiomeDASim]{mean_trend}} for details.
#' @param zero_trunc logical indicator designating whether simulated outcomes
#' should be zero truncated. default is set to TRUE
#' @param asynch_time logical indicator designed to randomly sample timepoints
#' over a specified interval if set to TRUE. default is FALSE.
#'
#' @importFrom graphics plot
#'
#' @examples
#' num_subjects_per_group <- 20
#' sim_obj <- mvrnorm_sim(n_control=num_subjects_per_group,
#'                        n_treat=num_subjects_per_group,
#'                        control_mean=5, sigma=1, num_timepoints=5,
#'                        t_interval=c(0, 4), rho=0.95, corr_str='ar1',
#'                        func_form='linear', beta=c(0, 0.25),
#'                        missing_pct=0.6, missing_per_subject=2)
#' #checking the output
#' head(sim_obj$df)
#'
#' #total number of observations is 2(num_subjects_per_group)(num_timeponts)
#' sim_obj$N
#'
#' #there should be approximately 60% of the IDs with missing observations
#' length(unique(sim_obj$miss_data$miss_id))/length(unique(sim_obj$df$ID))
#'
#' #checking the subject covariance structure
#' sim_obj$Sigma[seq_len(5), seq_len(5)]
#'
#' @return
#' This function returns a list with the following objects:
#'
#' \code{df} - data.frame object with complete outcome \code{Y}, subject ID,
#'     time, group, and outcome with missing data
#'
#' \code{Y} - vector of complete outcome
#'
#' \code{Mu} - vector of complete mean specifications used during simulation
#'
#' \code{Sigma} - block diagonal symmetric matrix of complete data used during
#'  simulation
#'
#' \code{N} - total number of observations
#'
#' \code{miss_data} - data.frame object that lists which ID's and timepoints
#' were randomly selected to induce missingness
#'
#' \code{Y_obs} - vector of outcome with induced missingness
#'
#' @export
mvrnorm_sim <- function(n_control, n_treat, control_mean, sigma, num_timepoints,
                        t_interval, rho, corr_str=c("ar1", "compound", "ind"),
                        func_form=c("linear", "quadratic", "cubic", "M", "W",
                                "L_up", "L_down"), beta, IP=NULL, missing_pct,
                        missing_per_subject, miss_val=NA, dis_plot=FALSE,
                        plot_trend=FALSE, zero_trunc=TRUE, asynch_time=FALSE) {
    corr_str <- match.arg(corr_str)
    func_form <- match.arg(func_form)
    n <- sum(n_control, n_treat)
    timepoints <- timepoint_process(num_timepoints, t_interval, n, asynch_time,
                                    missing_per_subject)
    N_c <- sum(timepoints$num_timepoints[seq_len(n_control)])
    t_tx <- timepoints$t[-seq_len(N_c)]
    trt_mean <- mean_trend(timepoints=t_tx, form=func_form,
                            beta=beta, IP=IP, plot_trend=plot_trend)
    treat_mean_mu <- control_mean + trt_mean$trend$mu
    mu_tot <- c(rep(control_mean, N_c), treat_mean_mu)
    rand_dt <- mvrnorm_corr_gen(n=n, obs=timepoints$num_timepoints,
                                t=timepoints$t, mu=mu_tot, sigma=sigma, rho=rho,
                                corr_str=corr_str, zero_trunc=zero_trunc)
    rand_dt$df$group <- ifelse(rand_dt$df$ID <= n_control, "Control",
                                "Treatment")
    missing_ids <- sample(seq_len(n), ceiling(n * missing_pct))
    missing_times <- unlist(lapply(missing_ids, function(x) {
        sample(2:num_timepoints, missing_per_subject)
    }))
    miss_data <- data.frame(miss_id=rep(missing_ids, each=missing_per_subject),
                            miss_time=missing_times)
    rand_dt$miss_data <- miss_data
    rand_dt$df$Y_obs <- rand_dt$df$Y
    rand_dt$Y_obs <- rand_dt$Y
    for (i in seq_len(nrow(miss_data))) {
        missing_ob <- which(rand_dt$df[, "ID"] == miss_data$miss_id[i] &
                                rand_dt$df[, "time"] == miss_data$miss_time[i])
        rand_dt$df[missing_ob, "Y_obs"] <- miss_val
        rand_dt$Y_obs[missing_ob] <- miss_val
    }
    if (dis_plot == TRUE) {
        p <- suppressWarnings(ggplot_spaghetti(y=rand_dt$df$Y_obs,
                            id=rand_dt$df$ID, time=rand_dt$df$time,
                            jit=0.1, group=rand_dt$df$group)) +
            labs(title="Simulated Microbiome Data from Multivariate Normal",
                    y="Normalized Reads", x="Time") +
            scale_linetype_manual(values=c("solid","dashed"), name="Group") +
            scale_color_manual(values=c("#F8766D", "#00BFC4"), name="Group")
        plot(p)
    }
    return(rand_dt)
}

#' Function for processing and checking the inputed timepoints
#'
#' To allow for increased flexibility the user may specify the number of
#' timepoints as either a single value or separately for each individual. There
#' is also an added option about whether to draw the timepoints evenly spaced
#' across the interval of interest or whether to randomly draw them.
#'
#' @param num_timepoints either an integer value specifying the number of
#'  timepoints per subject or a vector of timepoints for each subject. If
#'  supplying a vector the lenght of the vector must equal the total number of
#'  subjects.
#' @param t_interval numeric vector of length two specifying the interval of
#' time from which to draw observatoins \[t_1, t_q\]. Assumed to be equally
#' spaced over the interval unless \code{asynch_time} is set to TRUE.
#' @param n numeric value representing the total number of obserations
#' @param asynch_time logical indicator designed to randomly sample timepoints
#' over a specified interval if set to TRUE.
#'
#' @details
#' It is assummed that there is a known time interval of interest over which
#' samples will be collected longitudinally on subjects. This interval is
#' specified as \[t_1, t_q\]. All subjects are assumed to have baseline
#' observations, i.e., t_1.
#'
#' Over this study interval each subject can have a potentially different number
#' of measurements taken. In the most simple case we assume that all subjects
#' will have the same number of measurements and can specify
#' \code{num_timepoints} as a single scalar value. Otherwise, we must specify
#' how many timepoints will be collected for each individual. In this latter
#' case \code{num_timepoints} must have the same length as the number of
#' subjects.
#'
#' Finally, we can select whether we want the timepoints to be drawn at equal
#' spaces over our study interal, or whether we want to randomly sample
#' asynchronous timepoints. In the asynchronous case we randomly draw from a
#' uniform distribution over the study interval with the restriction that the
#' first observation must occur at t_1.
#'
#' @return
#' Returns a list of the number of timepoints and the times for each unit
#'
#' @keywords internal
timepoint_process <- function(num_timepoints, t_interval, n, asynch_time,
                              missing_per_subject){
    time_len <- length(num_timepoints)
    if(time_len != 1 && time_len != n){
        stop("length of timepoints misspecified", call.=FALSE)
    }
    if(time_len == 1){
        num_timepoints <- rep(num_timepoints, n)
    }
     if(!any(unlist(lapply(c(num_timepoints, t_interval), is.numeric)))){
        stop("num_timepoints and/or t_interval not numeric", call.=FALSE)
    }
    if(any(is.infinite(num_timepoints)) || any(num_timepoints<=0)){
        stop("num_timepoints has infinite or non positive values", call.=FALSE)
    }
    if(any(missing_per_subject > (num_timepoints - 1))){
        stop("Missing per subject greater than at least number of timepoints.",
             call.=FALSE)
    }
    if(length(t_interval)!=2){
        stop("time interval must be numeric vector of length 2", call.=FALSE)
    }
    if(any(is.infinite(t_interval))){
        stop("time interval must be finite", call.=FALSE)
    }
    if(asynch_time==FALSE){
        t <- lapply(num_timepoints, function(nt){
            seq(from=t_interval[1], to=t_interval[2], length.out=nt)
        })
    }else{
        t <- lapply(num_timepoints, function(nt){
            t <- runif(nt-1, min=t_interval[1], max=t_interval[2])
            t <- t[order(t)]
            t <- c(0, t)
        })
    }
    t <- unlist(t)
    return(list(num_timepoints=num_timepoints, t=t))
}

#' Generate Multivariate Random Normal Longitudinal Data
#'
#' For this methodology we assume that we draw a set of `n` independent each
#'  with \eqn{q_{i}} observations.
#'
#' @param n integer scalar representing the total number of individuals
#' @param obs integer or vector specifying the number of observations per
#'  indivdiual. If an integer then all indivdiuals are assummed to have the same
#'   number of obsevations. If a vector, then the vector must have length equal
#'    to \code{n} where each element specifies the number of observations for
#'     the \eqn{i^{th}} individual.
#' @param t vector corresponding to the timepoints for each individual.
#' @param mu integer or vector specifying the mean value for individuals.
#' If an integer then all indivdiuals are assummed to have the same mean.
#' If a vector, then the vector must have length equal to \code{n} where each
#'  element specifies the mean for the \eqn{i^{th}} individual.
#' @param sigma numeric scalar or vector specifying the standard deviation for
#'  observations.
#' @param rho numeric scalar value between \[0, 1\] specifying the amount of
#'  correlation between. assumes that the correlation is consistent for all
#'  subjects.
#' @param corr_str character value specifying the correlation structure.
#' Currently available methods are \'ar1\', \'compound\', and \'ind\' which
#' correspond to first-order autoregressive, compound or equicorrelation,
#'  and independence respecitvely.
#' @param zero_trunc logical value to specifying whether the generating
#'  distribution should come from a multivariate zero truncated normal or an
#'  untruncated multivariate normal. by default we assume that zero truncation
#'   occurs since this is assummed in our microbiome setting.
#'
#' @importFrom Matrix bdiag
#'
#' @examples
#' mvrnorm_corr_gen(n=15, obs=4, mu=20, sigma=2, rho=0.9, corr_str="ar1")
#'
#' @return
#' This function returns a list with the following objects:
#'
#' \code{df} - data.frame object with complete outcome \code{Y}, subject ID,
#'  time, group, and outcome with missing data
#'
#' \code{Y} - vector of complete outcome
#'
#' \code{Mu} - vector of complete mean specifications used during simulation
#'
#' \code{Sigma} - block diagonal symmetric matrix of complete data used during
#'  simulation
#'
#' \code{N} - total number of observations
#'
#' @export
mvrnorm_corr_gen <- function(n, obs, t, mu, sigma, rho,
                                corr_str=c("ar1", "compound", "ind"),
                                zero_trunc=TRUE) {
    if (!is.numeric(rho))
        stop("rho must be numeric", call.=FALSE)
    if (all.equal(n, as.integer(n)) != TRUE) {
        stop("n must be specified as an integer", call.=FALSE)
    }
    if (n <= 0)
        stop("n must be positive", call.=FALSE)
    if ((rho > 1 || rho < 0) && corr_str %in% c("ar1", "compound")) {
        stop("rho must be in [0, 1]", call.=FALSE)
    }
    if (corr_str == "ind" && is.null(rho) == FALSE) {
        warning("ignoring rho coefficient", call.=FALSE)
    }
    N <- sum(obs)
    s_df <- data.frame(time=t, id=rep(seq_len(n), obs))
    id_split <- split(s_df, s_df$id)
    block_diag_list <- lapply(id_split, function(s){
        sigma_corr_function(t=s$time, sigma, corr_str,rho)
        })
    Sigma <- bdiag(block_diag_list)
    Y <- trunc_bugs(Y, N, Mu=mu, Sigma, zero_trunc)
    df <- data.frame(Y, ID=s_df$id, time=s_df$time)
    return(list(df=df, Y=Y, Mu=mu, Sigma=Sigma, N=N))
}

#' Checking input for scalar or vector valued
#'
#' This function allows users to specify just a scalar value as a parameter in
#'  the generating function, and here it is converted to a vector
#' of the proper length.
#'
#' Note that `n` refers to the number of individuals in the
#' \code{\link[microbiomeDASim]{mvrnorm_corr_gen}}
#'
#' @param input a variable that will be supplied that we want to either check
#'  the length or replicate to the specified length
#' @param n an integer value specifying the desired vector length if a scalar
#'  is provided
#' @param N an integer value specifying the total number of observations across
#' individuals.
#'
#' @keywords internal
#'
#' @return
#' return a vector that has same length as specified n
vector_scalar_check <- function(input, n, N) {
    if (length(input) == 1) {
        output <- rep(input, N)
    } else if (length(input) == n | length(input) == N) {
        output <- input
    } else{
        stop("incorrect vector length", call.=FALSE)
    }
    return(output)
}

#' Function for inducing truncation of outcome
#'
#' @param Y The original N x 1 vector of simulated multivariate outcomes
#' @param N Total number of observations equal to sum of repeated measurements
#' for all individuals
#' @param Mu N x 1 vector representing the mean values
#' @param Sigma N x N numeric matrix representing the covariance matrix for
#' the feature
#' @param zero_trunc Logical indicator whether to perform zero-truncation
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom MASS mvrnorm
#'
#' @keywords internal
#'
#' @return
#' Potentially truncated outcome vector Y
trunc_bugs <- function(Y, N, Mu, Sigma, zero_trunc){
    if (zero_trunc == TRUE) {
        if (N < 1000) {
            accep_prob <- pmvnorm(lower=rep(0, N), mean=Mu,
                                            sigma=as.matrix(Sigma))
            if (accep_prob[1] > 0.1) {
                Y <- as.numeric(rtmvnorm(1, mean=Mu, sigma=Sigma,
                                        lower=rep(0, N), algorithm="rejection"))
            } else {
                Y <- mvrnorm(1, Mu, Sigma)
                Y <- ifelse(Y < 0, 0, Y)
            }
        } else {
            Y <- mvrnorm(1, Mu, Sigma)
            Y <- ifelse(Y < 0, 0, Y)
        }
    } else {
        Y <- mvrnorm(1, Mu, Sigma)
    }
}

#' Generating the longitudinal correlation matrix for repeated observations
#'
#' @param t timepoints for repeated observations
#' @param sigma the standard deviation parameter for the covariance matrix
#' @param corr_str the type of correlatin structure chosen. options currently
#' available include "ar1", "compound", and "ind"
#' @param rho the correlation coefficient for non-independent structures
#'
#' @keywords internal
#'
#' @return
#' Return the covariance matrix V as a list
sigma_corr_function <- function(t, sigma, corr_str, rho){
    if (corr_str == "ar1") {
        H <- abs(outer(t, t, "-"))
        V <- sigma * rho^H
        p <- nrow(V)
        pl <- seq_len(p)
        V[cbind(pl, pl)] <- V[cbind(pl, pl)] * sigma
    } else if (corr_str == "compound") {
        H <- outer(rep(1, length(t)), rep(1, length(t))) - diag(1, length(t))
        V <- sigma * rho^H
        p <- nrow(V)
        pl <- seq_len(p)
        V[cbind(pl, pl)] <- V[cbind(pl, pl)] * sigma
    } else if (corr_str == "ind") {
        V <- sigma * diag(1, length(t))
        p <- nrow(V)
        pl <- seq_len(p)
        V[cbind(pl, pl)] <- V[cbind(pl, pl)] * sigma
    }
    return(V)
}
