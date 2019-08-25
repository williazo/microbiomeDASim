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
#' @param num_timepoints integer value specifying the number of timepoints per
#'  subject.
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
#'
#' @importFrom graphics plot
#'
#' @examples
#' num_subjects_per_group <- 20
#' sim_obj <- mvrnorm_sim(n_control=num_subjects_per_group,
#'                        n_treat=num_subjects_per_group,
#'                        control_mean=5, sigma=1, num_timepoints=5,
#'                        rho=0.95, corr_str='ar1', func_form='linear',
#'                        beta=c(0, 0.25),
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
                        rho, corr_str=c("ar1", "compound", "ind"),
                        func_form=c("linear", "quadratic", "cubic", "M", "W",
                                "L_up", "L_down"), beta, IP=NULL, missing_pct,
                        missing_per_subject, miss_val=NA, dis_plot=FALSE,
                        plot_trend=FALSE, zero_trunc=TRUE) {
    corr_str <- match.arg(corr_str)
    func_form <- match.arg(func_form)
    if (missing_per_subject > (num_timepoints - 1)) {
        stop("Value of missing_per_subject > than num_timepoints-1.",
                call.=FALSE)
    }
    n <- sum(n_control, n_treat)
    trt_mean <- mean_trend(timepoints=seq_len(num_timepoints), form=func_form,
                            beta=beta, IP=IP, plot_trend=plot_trend)
    treat_mean_mu <- control_mean + trt_mean$trend$mu
    mu_tot <- c(rep(control_mean, n_control * num_timepoints),
                rep(treat_mean_mu, n_treat))
    rand_dt <- mvrnorm_corr_gen(n=n, obs=num_timepoints, mu=mu_tot, sigma=sigma,
                                rho=rho, corr_str=corr_str,
                                zero_trunc=zero_trunc)
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
mvrnorm_corr_gen <- function(n, obs, mu, sigma, rho,
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
    obs <- vector_scalar_check(obs, n)
    mu <- vector_scalar_check(mu, n)
    sigma <- vector_scalar_check(sigma, n)
    corr_str <- vector_scalar_check(corr_str, n)
    rho <- vector_scalar_check(rho, n)
    N <- sum(obs)
    if (length(mu) == n) {
        Mu <- unlist(mapply(x=mu, y=obs, function(x, y) {
            mu_i <- rep(x, y)
            return(list(mu_i))
        }))
    } else if (length(mu) == N) {
        Mu <- mu
    } else {
        stop("mu must be a scalar, vector of length n, or vector of length N",
                call.=FALSE)
    }
    block_diag_list <- mapply(obs, sigma, corr_str, rho,FUN=sigma_corr_function)
    Sigma <- bdiag(block_diag_list)
    Y <- trunc_bugs(Y, N, Mu, Sigma, zero_trunc)
    df <- data.frame(Y, ID=rep(seq_len(n), obs),
                        time=unlist(lapply(obs, seq_len)))
    return(list(df=df, Y=Y, Mu=Mu, Sigma=Sigma, N=N))
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
#'
#' @return
#' return a vector that has same length as specified n
vector_scalar_check <- function(input, n) {
    if (length(input) == 1) {
        output <- rep(input, n)
    } else if (length(input) > 1) {
        output <- input
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
#' @param obs number of repeated observations per individual
#' @param sigma the standard deviation parameter for the covariance matrix
#' @param corr_str the type of correlatin structure chosen. options currently
#' available include "ar1", "compound", and "ind"
#' @param rho the correlation coefficient for non-independent structures
#'
#' @return
#' Return the covariance matrix V as a list
sigma_corr_function <- function(obs, sigma, corr_str, rho){
    if (corr_str == "ar1") {
        H <- abs(outer(seq_len(obs), seq_len(obs), "-"))
        V <- sigma * rho^H
        p <- nrow(V)
        pl <- seq_len(p)
        V[cbind(pl, pl)] <- V[cbind(pl, pl)] * sigma
    } else if (corr_str == "compound") {
        H <- outer(rep(1, obs), rep(1, obs)) - diag(1, obs)
        V <- sigma * rho^H
        p <- nrow(V)
        pl <- seq_len(p)
        V[cbind(pl, pl)] <- V[cbind(pl, pl)] * sigma
    } else if (corr_str == "ind") {
        V <- sigma * diag(1, obs)
        p <- nrow(V)
        pl <- seq_len(p)
        V[cbind(pl, pl)] <- V[cbind(pl, pl)] * sigma
    }
    return(list(V))
}
