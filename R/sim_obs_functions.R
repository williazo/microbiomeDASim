#' Simulate Microbiome Longitudinal Data from Multivariate Random Normal with
#' Observed Data
#'
#' This function is used in the
#' \code{\link[microbiomeDASim]{gen_norm_microbiome_obs}} call.
#'
#' @param id vector of length \code{N} that identifies repeated measurements for
#' each unit
#' @param time vector of length \code{N} that determines when values will be
#' sampled for each unit
#' @param group factor vector with two levels indicating the group assignment
#' for each respective id
#' @param ref character value identifying which group value to treat as control
#' and which value to treat as treatment
#' @param control_mean numeric value specifying the mean value for control
#' subjects. all control subjects are assummed to have the same population mean
#'  value.
#' @param sigma numeric value specifying the global population standard
#'  deviation for both control and treated individuals.
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
#' @param dis_plot logical argument on whether to plot the simulated data or
#'     not. by default plotting is turned off.
#' @param plot_trend specifies whether to plot the true mean trend. see
#' \code{\link[microbiomeDASim]{mean_trend}} for details.
#' @param zero_trunc logical indicator designating whether simulated outcomes
#' should be zero truncated. default is set to TRUE
#'
#' @importFrom graphics plot
#' @importFrom stats relevel
#'
#' @examples
#' set.seed(011520)
#' id_list <- lapply(seq_len(30), function(i){
#' obs <- sample(seq_len(10), size=1)
#' id_rep <- rep(i, obs)
#' })
#'
#' time_interval <- c(0, 10)
#' time_list <- lapply(id_list, function(x){
#' time_len <- length(x)
#' times <- runif(time_len, min=time_interval[1], max=time_interval[2])
#' times <- times[order(times)]
#' })
#'
#' group_list <- lapply(id_list, function(x){
#' group_len <- length(x)
#' tx_ind <- sample(seq_len(2), 1)
#' tx_group <- ifelse(tx_ind==1, "Control", "Treatment")
#' groups <- rep(tx_group, group_len)
#' })
#' id <- unlist(id_list)
#' group <- factor(unlist(group_list), levels = c("Control", "Treatment"))
#' time <- unlist(time_list)
#'
#' # N=173 total repeated measurements
#' length(id)
#'
#' # 15 control and 15 treated subjects
#' table(group[unique(id)])
#'
#' # control times
#' ct <- unlist(lapply(unique(id[group=="Control"]), function(x){
#' length(id[id==x])
#' }))
#'
#' #treatment times
#' tt <- unlist(lapply(unique(id[group=="Treatment"]), function(x){
#' length(id[id==x])
#' }))
#'
#' # on average the treatment group has one more observation than control
#' mean(ct)
#' mean(tt)
#'
#'mvrnorm_sim_obs(id=id, time=time, group=group, ref="Control", control_mean=2,
#'                sigma=1, rho=0.7, corr_str="compound", func_form="L_up",
#'                beta=1, IP=5, plot_trend=TRUE, dis_plot=TRUE, zero_trunc=TRUE)
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
#' @export
mvrnorm_sim_obs <- function(id, time, group, ref, control_mean, sigma, rho,
        corr_str=c("ar1", "compound", "ind"),
        func_form=c("linear", "quadratic", "cubic", "M", "W", "L_up", "L_down"),
        beta, IP=NULL, dis_plot=FALSE, plot_trend=FALSE, zero_trunc=TRUE){
    corr_str <- match.arg(corr_str)
    func_form <- match.arg(func_form)
    if(!is.factor(group)) stop("group must be factor variable", call.=FALSE)
    if(length(levels(group))!=2) stop("group must be binary", call.=FALSE)
    if(!ref%in%levels(group)){
        stop("ref must be specified as one of group values", call.=FALSE)}
    if(length(id)!=length(time) || length(id)!=length(group)
        || length(time)!=length(group)){
        stop("id, time, and group must all be length N", call.=FALSE)
    }
    N <- length(id)
    n <- length(unique(id))
    df <- data.frame(id, time, group)
    df$group <- relevel(df$group, ref=ref)
    df_group <- split(df, group)
    trt_mean <- mean_trend(timepoints=df_group[[2]]$time, form=func_form,
                            beta=beta, IP=IP, plot_trend=plot_trend)
    df_group[[2]]$mu <- control_mean + trt_mean$trend$mu
    df_group[[1]]$mu <- control_mean
    df <- data.frame(do.call(rbind, df_group))
    ob_ct <- unlist(lapply(unique(df$id), function(x) length(df$id[df$id==x])))
    rand_dt <- mvrnorm_corr_gen(n=n, obs=ob_ct, t=df$time, mu=df$mu,
                                sigma=sigma, rho=rho, corr_str=corr_str,
                                zero_trunc=zero_trunc)
    rand_dt$df$ID <- df$id
    rand_dt$df$group <- df$group
    if (dis_plot == TRUE) {
        p <- suppressWarnings(ggplot_spaghetti(y=rand_dt$df$Y,id=rand_dt$df$ID,
                    time=rand_dt$df$time, jit=0.1, group=rand_dt$df$group)) +
            labs(title="Simulated Microbiome Data from Multivariate Normal",
                    y="Normalized Reads", x="Time") +
            scale_linetype_manual(values=c("solid","dashed"), name="Group") +
            scale_color_manual(values=c("#F8766D", "#00BFC4"), name="Group")
        plot(p)
    }
    return(rand_dt)
}


#' Generate Longitduinal Differential Abundance from Multivariate Normal with
#' Observed Data
#'
#' @param features numeric value specifying the number of features/microbes to
#'     simulate. Default is 10.
#' @param diff_abun_features numeric value specifying the number of
#'     differentially abundant features. Default is 5.
#' @param id vector of length \code{N} that identifies repeated measurements for
#' each unit
#' @param time vector of length \code{N} that determines when values will be
#' sampled for each unit
#' @param group factor vector with two levels indicating the group assignment
#' for each respective id
#' @param ref character value identifying which group value to treat as control
#' and which value to treat as treatment
#' @param sigma numeric value specifying the global population standard
#'  deviation for both control and treated individuals.
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
#' @param dis_plot logical argument on whether to plot the simulated data or
#'     not. by default plotting is turned off.
#' @param plot_trend specifies whether to plot the true mean trend. see
#' \code{\link[microbiomeDASim]{mean_trend}} for details.
#' @param zero_trunc logical indicator designating whether simulated outcomes
#' should be zero truncated. default is set to TRUE
#'
#' @importFrom pbapply pblapply
#'
#' @return This function returns a list with the following objects
#'
#' \code{Y} The full simulated feature sample matrix where each row represent a
#'     feature and each column a sample. Note that the differential and
#'     non-differential bugs are marked by row.names
#'
#' @examples
#' set.seed(011520)
#' id_list <- lapply(seq_len(60), function(i){
#' obs <- sample(5:10, size=1)
#' id_rep <- rep(i, obs)
#' })
#'
#' time_interval <- c(0, 10)
#' time_list <- lapply(id_list, function(x){
#' time_len <- length(x)
#' times <- runif(time_len, min=time_interval[1], max=time_interval[2])
#' times <- times[order(times)]
#' })
#'
#' group_list <- lapply(id_list, function(x){
#' group_len <- length(x)
#' tx_ind <- sample(seq_len(2), 1)
#' tx_group <- ifelse(tx_ind==1, "Control", "Treatment")
#' groups <- rep(tx_group, group_len)
#' })
#' id <- unlist(id_list)
#' group <- factor(unlist(group_list), levels = c("Control", "Treatment"))
#' time <- unlist(time_list)
#'
#' # control times
#' ct <- unlist(lapply(unique(id[group=="Control"]), function(x){
#' length(id[id==x])
#' }))
#'
#' tt <- unlist(lapply(unique(id[group=="Treatment"]), function(x){
#' length(id[id==x])
#' }))
#' mean(ct)
#' mean(tt)
#'
#' gen_norm_microbiome_obs(features=4, diff_abun_features=2,
#' id=id, time=time, group=group, ref="Control", control_mean=2,
#'                sigma=1, rho=0.7, corr_str="compound", func_form="L_up",
#'                beta=1, IP=5, zero_trunc=TRUE)
#'
#' @export
gen_norm_microbiome_obs <- function(features=10, diff_abun_features=5, id, time,
        group, ref, control_mean, sigma, rho,
        corr_str=c("ar1", "compound", "ind"),
        func_form=c("linear", "quadratic", "cubic", "M", "W", "L_up", "L_down"),
        beta, IP=NULL, dis_plot=FALSE, plot_trend=FALSE, zero_trunc=TRUE){
    gen_microbiome_norm_feature_check(features, diff_abun_features)
    no_diff_feat <- features - diff_abun_features
    if(diff_abun_features > 0){
        message("Simulating Diff Bugs\n")
        diff_bugs <- pblapply(seq_len(diff_abun_features), function(x){
            mvrnorm_sim_obs(id, time, group, ref, control_mean, sigma, rho,
                            corr_str, func_form, beta, IP, dis_plot,
                            plot_trend, zero_trunc)})
        diff_Y <- lapply(diff_bugs, function(x) return(x$Y))
        diff_Y <- matrix(unlist(diff_Y), nrow=diff_abun_features, byrow=TRUE)
    } else{
        diff_Y <- NULL
        diff_bugs <- NULL
    }
    if (no_diff_feat > 0) {
        message("Simulating No-Diff Bugs\n")
        nodiff_bugs <- pblapply(seq_len(no_diff_feat), function(x){
            mvrnorm_sim_obs(id, time, group, ref, control_mean, sigma, rho,
                        corr_str, func_form="linear", beta=c(0, 0), IP,
                        dis_plot, plot_trend, zero_trunc)})
        null_Y <- lapply(nodiff_bugs, function(x) return(x$Y))
        null_Y <- matrix(unlist(null_Y), nrow=no_diff_feat, byrow=TRUE)
    } else{
        null_Y <- NULL
        nodiff_bugs <- NULL
    }
    final_output <- final_output_gen(no_diff_feat, diff_abun_features,
                                        diff_Y, null_Y, diff_bugs, nodiff_bugs,
                                        final_output=NULL)
    return(final_output)
}
