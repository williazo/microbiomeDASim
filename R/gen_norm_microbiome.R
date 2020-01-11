#' Generate Longitduinal Differential Abundance from Multivariate Normal
#'
#' @param features numeric value specifying the number of features/microbes to
#'     simulate. Default is 10.
#' @param diff_abun_features numeric value specifying the number of
#'     differentially abundant features. Default is 5.
#' @param n_control integer value specifying the number of control individuals
#' @param n_treat integer value specifying the number of treated individuals
#' @param control_mean numeric value specifying the mean value for control
#' subjects. all control subjects are assummed to have the same population mean
#'  value.
#' @param sigma numeric value specifying the global population standard
#'  deviation for both control and treated individuals.
#' @param num_timepoints integer value specifying the number of timepoints per
#'  subject.
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
#' @importFrom pbapply pblapply
#'
#' @return This function returns a list with the following objects
#'
#' \code{Y} The full simulated feature sample matrix where each row represent a
#'     feature and each column a sample. Note that the differential and
#'     non-differential bugs are marked by row.names
#'
#' @examples
#'gen_norm_microbiome(features = 5, diff_abun_features = 2,
#'                n_control = 10, n_treat = 10, control_mean = 8, sigma = 1,
#'                num_timepoints = 5, t_inverval=c(0, 4), rho = 0.8,
#'                corr_str = "compound", func_form = "linear", beta =  c(0, 1),
#'                missing_pct = 0.3, missing_per_subject = 2)
#'
#' @export
gen_norm_microbiome <- function(features=10, diff_abun_features=5,
                                n_control, n_treat, control_mean,
                                sigma, num_timepoints, t_interval, rho,
                                corr_str=c("ar1", "compound", "ind"),
                                func_form=c("linear", "quadratic", "cubic",
                                            "M", "W", "L_up", "L_down"),
                                beta, IP=NULL, missing_pct, missing_per_subject,
                                miss_val=NA, dis_plot=FALSE, plot_trend=FALSE,
                                zero_trunc=TRUE, asynch_time=FALSE) {
    gen_microbiome_norm_feature_check(features, diff_abun_features)
    no_diff_feat <- features - diff_abun_features
    if(diff_abun_features > 0){
        message("Simulating Diff Bugs\n")
        diff_bugs <- pblapply(seq_len(diff_abun_features), function(x){
            mvrnorm_sim(n_control, n_treat, control_mean, sigma, num_timepoints,
                        t_interval, rho, corr_str, func_form, beta, IP,
                        missing_pct, missing_per_subject, miss_val, dis_plot,
                        plot_trend, zero_trunc, asynch_time)})
        diff_Y <- lapply(diff_bugs, function(x) return(x$Y))
        diff_Y <- matrix(unlist(diff_Y), nrow=diff_abun_features, byrow=TRUE)
    } else{
        diff_Y <- NULL
        diff_bugs <- NULL
    }
    if (no_diff_feat > 0) {
        message("Simulating No-Diff Bugs\n")
        nodiff_bugs <- pblapply(seq_len(no_diff_feat), function(x){
            mvrnorm_sim(n_control, n_treat, control_mean, sigma, num_timepoints,
                        t_interval, rho, corr_str, func_form, beta, IP,
                        missing_pct, missing_per_subject, miss_val, dis_plot,
                        plot_trend, zero_trunc, asynch_time)})
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

#' Checking that features are specified appopriately
#'
#' @param features Numeric value specifying the total number of features to
#' simulate in the microbiome. Must be greater than zero
#' @param diff_abun_features Number of features to simulate with differentially
#' abundant pattern. Must be between zero and number of features specified
#'
#' @keywords internal
#'
#' @return
#' Potential warning message if no differentially abundant features or all
#' differentially abundant features are specified
gen_microbiome_norm_feature_check <- function(features, diff_abun_features){
    if(features <= 0 ){
        stop("must specify features greater than zero", call.=FALSE)
    }
    if (diff_abun_features > features) {
        stop("differential abundance must be <= to total features",
                call.=FALSE)
    } else if (diff_abun_features == features) {
        warning("all features will be simulated with differential abundance",
                call.=FALSE)
    } else if (diff_abun_features == 0) {
        warning("no differential abundance features specified",
                call.=FALSE)
        warning("ignoring functional form and beta specification",
                call.=FALSE)
    }
}

#' Generating the final combined bug output
#'
#' @param no_diff_feat number of non differentially abundant features
#' @param diff_abun_features number of differentially abundant features
#' @param diff_Y simulated outcome for differentially abundant features
#' @param null_Y simulated outcome for non differentially abundant features
#' @param diff_bugs sample information for differentially abundant features
#' @param nodiff_bugs sample information for non differentially abundant
#' features
#' @param final_output final object that will store the simulated data
#'
#' @keywords internal
#'
#' @return
#' final output list with the OTU table and corresponding bug feature data.frame
final_output_gen <- function(no_diff_feat, diff_abun_features, diff_Y,
                            null_Y, diff_bugs, nodiff_bugs, final_output=NULL){
    if(no_diff_feat > 0  && diff_abun_features > 0){
        final_output$Y <- rbind(diff_Y, null_Y)
        s_id <- paste0("Sample_", rownames(diff_bugs[[1]]$df))
        colnames(final_output$Y) <- s_id
        rownames(final_output$Y) <- c(paste0("Diff_Bug",
                                                seq_len(diff_abun_features)),
                                        paste0("NoDiffBug_",
                                                seq_len(no_diff_feat)))
        Y_vars <- grep("Y", names(diff_bugs[[1]]$df))
        final_output$bug_feat <- diff_bugs[[1]]$df[, -Y_vars]
        final_output$bug_feat$Sample_ID <- s_id
        }
    else if(no_diff_feat == 0 && diff_abun_features > 0){
        final_output$Y <- diff_Y
        s_id <- paste0("Sample_", rownames(diff_bugs[[1]]$df))
        colnames(final_output$Y) <- s_id
        rownames(final_output$Y) <- paste0("Diff_Bug",
                                            seq_len(diff_abun_features))
        Y_vars <- grep("Y", names(diff_bugs[[1]]$df))
        final_output$bug_feat <- diff_bugs[[1]]$df[, -Y_vars]
        final_output$bug_feat$Sample_ID <- s_id
        }
    else{
        final_output$Y <- null_Y
        s_id <- paste0("Sample_", rownames(nodiff_bugs[[1]]$df))
        colnames(final_output$Y) <- s_id
        rownames(final_output$Y) <- paste0("NoDiff_Bug", seq_len(no_diff_feat))
        Y_vars <- grep("Y", names(nodiff_bugs[[1]]$df))
        final_output$bug_feat <- nodiff_bugs[[1]]$df[, -Y_vars]
        final_output$bug_feat$Sample_ID <- s_id
        }
    return(final_output)
}
