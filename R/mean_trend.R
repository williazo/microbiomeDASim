#' Function for Generating Various Longitudinal Mean Trends
#'
#' In order to investigate different functional forms of longitudinal
#' differential abundance we allow the mean time trend to take a variety of
#'  forms.
#' These functional forms include linear, quadratic, cubic, M, W, L_up, or
#'  L_down. For each form the direction/concavity/fold change can be specified
#'   using the beta parameter.
#' @param timepoints numeric vector specifying the points to fit the functional
#'  trend.
#' @param form character value specifying the type of time trend. Options
#'  include 'linear', 'quadratic', 'cubic', 'M', 'W', 'L_up', and 'L_down'.
#' @param beta vector specifying the appropriate parameters for the equation.
#' In the case of 'linear', beta should be a two-dimensional vector specifying
#' the intercept and slope. See details for the further explanation of the beta
#' value for each form.
#' @param IP vector specifying the inflection points where changes occur for
#' functional forms M, W, and L trends.
#' @param plot_trend logical value indicating whether a plot should be produced
#' for the time trend. By default this is set to TRUE.
#'
#'
#' @details
#' Linear Form Notes:
#' \deqn{f(x)=\beta_0+\beta_{1}x+\beta_{2}x^2}
#' \itemize{
#'  \item{Sign of \eqn{\beta_1} determines whether the trend is increasing (+)
#'  or decreasing (-)}
#'  }
#'
#' Quadratic Form Notes:
#' \deqn{f(x)=\beta_0+\beta_{1}x+\beta_{2}x^2}
#' \itemize{
#'  \item{Critical point for quadratic function occurs at the point
#'   \eqn{\frac{-\beta_1}{2\beta_2}}}
#'  \item{\eqn{\beta_2} determines whether the quadratic is concave up (+) or
#'   concave down (-)}
#'  }
#'
#' Cubic Form Notes:
#' \deqn{f(x)=\beta_0+\beta_{1}x+\beta_{2}x^2+\beta_{3}x^{3}}
#' \itemize{
#'  \item{Point of Inflection for cubic function occurs
#'   \eqn{\frac{-\beta_{2}}{(3\beta_{3})}}}
#'  \item{Critical points for cubic function occur at
#'   \eqn{\frac{-\beta_{2}\pm\sqrt{\beta_2^{2}-3\beta_{1}\beta_{3}}}{3\beta_3}}}
#'  \item{Can generate piecewise linear trends, i.e. 'V' form, by placing either
#'   one of the IP points outside of the timepoints specified}
#' }
#'
#' M/W Form Notes:
#' \itemize{
#'  \item{Must specify beta as (\eqn{\beta_0}, \eqn{\beta_1}) and IP
#'  as (\eqn{IP_1}, \eqn{IP_2}, \eqn{IP_3})}
#'  \item{This form should be specified with an initial intercept,
#'  \eqn{\beta_0}, and slope, \eqn{\beta_1},
#'  that will connect to the first point of change (IP) specified.}
#'  \item{Subsequent slopes are constructed such that the mean value at the
#'   second IP value and final timepoint are 0}
#'  \item{The mean value at the third IP is set to be equal to the calculcated
#'   mean value at the first IP based on the specified intercept and slope.}
#'  \item{\eqn{\beta_0}=intercept, i.e. timepoint when y=0}
#'  \item{\eqn{\beta_1}=slope between \eqn{\beta_0} and \eqn{IP_1}}
#' }
#'
#' L_up Form Notes:
#'
#' The structure of this form assumes that there is no trend from \eqn{t_{1}} to
#'  \eqn{IP_{1}}.
#' Then at the point of change specified, \eqn{IP_{1}}, there occurs a linearly
#'  increasing trend with slope equal to \eqn{\beta_{slope}} up to the last
#'  specified timepoint \eqn{t_{q}}.
#' \itemize{
#'  \item{Must specify beta as (\eqn{\beta_{slope}}), and must be positive}
#'  \item{Specify a single point of change (IP) variable where positive trend
#'   will start}
#'  \item{IP must be between \[\eqn{t_{1}}, \eqn{t_{q}}\]}
#' }
#'
#' L_down Form Notes:
#'
#' Similarily, the L_down form assumes that there are two region within the
#'  range of timepoints. The first region is a decreasing trend and the second
#'   region has no trend.
#' The decreasing trend must start with a Y intercept greater than zero, and the
#'  slope must be specified as negative. There is one point of change (IP),
#'   but this is
#' calculated automatically based on the values of the Y intercept and slope
#' provided, IP=\eqn{-\beta_{yintercept}/\beta_{slope}}.
#' \itemize{
#'  \item{Must specify beta as (\eqn{\beta_{yintercept}}, \eqn{\beta_{slope}})
#'   where \eqn{\beta_{yintercept}}>0 and \eqn{\beta_{slope}}<0}
#'  \item{IP variable should be specified as NULL, if value is provided it will
#'   be ignored.}
#' }
#'
#' @examples
#' #Quadratic Form
#' mean_trend(timepoints=seq(0, 6, length.out=20),
#'                form='quadratic', beta=1/4 * c(-1, 3, -0.5), plot_trend=TRUE)
#' #M Form
#' mean_trend(timepoints=seq(0, 10,length.out=100), form='M',
#'                beta=c(0, 5), IP=10 * c(1/4, 2/4, 3/4), plot_trend=TRUE)
#' #in this case the IP points are selected so that peaks are evenly
#' #distributed but this does not have to be true in general
#'
#' #L_up Form
#' mean_trend(timepoints=seq(0, 10, length.out=100), form='L_up',
#'            beta=1, IP=5, plot_trend=TRUE)
#'
#' #L_down Form
#' mean_trend(timepoints=seq(0, 10,length.out=100), form='L_down',
#'            beta=c(4, -0.5), IP=NULL, plot_trend=TRUE)
#'
#' @return
#' This function returns a list of the following
#'
#' \code{form} - character value repeating the form selected
#'
#' \code{trend} - data.frame with the variables \code{mu} representing the
#' estimated mean value at \code{timepoints} used for fitting the trend
#'
#' \code{beta} - returning the numeric vector used to fit the functional form
#'
#'
#' @import ggplot2
#'
#' @export
mean_trend <- function(timepoints, form=c("linear", "quadratic", "cubic",
                                            "M", "W", "L_up", "L_down"),
                        beta, IP=NULL, plot_trend=FALSE) {
    form <- match.arg(form)
    if (!is.logical(plot_trend)) {
        stop("plot_trend must be logical", call.=FALSE)
    }
    form_beta_check(form, beta, IP, timepoints)
    IP <- IP_form_check(form, beta, IP, timepoints)
    beta <- mean_trend_beta_vec(form, beta, IP, timepoints)
    design_mat <- mean_trend_design_mat(form, beta, IP, timepoints)
    mu <- design_mat %*% beta
    fit_df <- data.frame(mu=mu, timepoints=timepoints)
    if (plot_trend == TRUE) {
        print(ggplot(fit_df, aes(x=timepoints, y=mu)) +
                    geom_point() +
                    geom_line())
    }
    return(list(form=form, trend=fit_df, beta=beta))
}

#' Beta Specification Check
#'
#' Function for checking that the appopriate beta parameters are specified for
#' each of the mean trend specifications
#'
#' @param form character value specifying the type of time trend. Options
#'  include 'linear', 'quadratic', 'cubic', 'M', 'W', 'L_up', and 'L_down'.
#' @param beta vector specifying the appropriate parameters for functional
#' trend. See details of \code{\link{mean_trend}} for explanation for each form
#' @param IP vector specifying the inflection points. See details of
#' \code{\link{mean_trend}} for explanation for each form
#' @param timepoints numeric vector specifying the points to fit the functional
#'  trend.
#'
#' @return
#' Nothing returned unless an error is returned.
form_beta_check <- function(form, beta, IP, timepoints){
    form_lbl <- factor(form, levels=c("linear", "quadratic", "cubic",
                                        "M", "W", "L_up", "L_down"))
    expected_beta <- c(2, 3, 4, 2, 2, 1, 2)
    beta_list <- list(c("beta_0", "beta_1"),
                        c("beta_0", "beta_1", "beta_2"),
                        c("beta_0", "beta_1", "beta_2", "beta_3"),
                        c("beta_0", "beta_1"),
                        c("beta_0", "beta_1"),
                        c("beta_slope"),
                        c("beta_yintercept", "beta_slope"))
    if (length(beta) != expected_beta[as.numeric(form_lbl)]) {
        stop(paste0("For a ", form, " trend, must specify beta as (",
                    paste(beta_list[[as.numeric(form_lbl)]],
                            collapse=", "), ")."), call.=FALSE)
    }
    if (form == "M" || form == "W"){
        if (beta[2] > 0 && form == "W") {
            warning("Second beta term should be negative", call.=FALSE)
        }
        if (beta[2] < 0 && form == "M") {
            warning("Second beta term should be positive", call.=FALSE)
        }
    }
    else if (form == "L_down"){
        if (beta[2] >= 0) {
            stop("Second element of beta vector be negative", call.=FALSE)
        }
        if (beta[1] <= 0) {
            stop("First element of beta vector must be positive", call.=FALSE)
        }
    }

}



#' Create Design Matrix for \code{\link{mean_trend}} function
#'
#' By taking in the user specified parameters, we can return a design matrix
#' to use when creating the differential longitudinal abundance.
#' @param form character value specifying the type of time trend. Options
#'  include 'linear', 'quadratic', 'cubic', 'M', 'W', 'L_up', and 'L_down'.
#' @param beta vector specifying the appropriate parameters for functional
#' trend. See details of \code{\link{mean_trend}} for explanation for each form
#' @param IP vector specifying the inflection points. See details of
#' \code{\link{mean_trend}} for explanation for each form
#' @param timepoints numeric vector specifying the points to fit the functional
#'  trend.
#'
#' @return
#' Numeric matrix with values that will be used to generate functional trends

mean_trend_design_mat <- function(form, beta, IP, timepoints){
    if (form == "linear") {
        design_mat <- cbind(rep(1, length(timepoints)),
                            timepoints)
    } else if (form == "quadratic") {
        design_mat <- cbind(intercept=rep(1, length(timepoints)),
                            t=timepoints,
                            t_2=timepoints^2)
    } else if (form == "cubic") {
        design_mat <- cbind(intercept=rep(1, length(timepoints)),
                            t=timepoints,
                            t_2=timepoints^2,
                            t_3=timepoints^3)
    } else if (form == "M" || form == "W") {
        design_mat <- cbind(rep(1, length(timepoints)),
                            I(timepoints < IP[1]) * timepoints,
                            I(timepoints >= IP[1]
                                & timepoints < IP[2]) |
                                I(timepoints >= IP[3]),
                            I(timepoints >= IP[1]
                                & timepoints < IP[2]) * (timepoints - IP[1]),
                            I(timepoints >= IP[2]
                                & timepoints < IP[3]) * (timepoints - IP[2]),
                            I(timepoints >= IP[3]) * (timepoints - IP[3]))
    } else if (form == "L_up") {
        design_mat <- cbind(I(timepoints >= IP),
                            I(timepoints >= IP) * timepoints)
    } else if (form == "L_down") {
        design_mat <- cbind(I(timepoints < IP),
                            I(timepoints < IP) * timepoints)

    }
    return(design_mat)
}

#' Create beta vector for \code{mean_trend} for all functional forms
#'
#' @param form character value specifying the type of time trend. Options
#'  include 'linear', 'quadratic', 'cubic', 'M', 'W', 'L_up', and 'L_down'.
#' @param beta vector specifying the appropriate parameters for functional
#' trend. See details of \code{\link{mean_trend}} for explanation for each form
#' @param IP vector specifying the inflection points. See details of
#' \code{\link{mean_trend}} for explanation for each form
#' @param timepoints numeric vector specifying the points to fit the functional
#'  trend.
#'
#' @return
#' Vector with beta values used to create mean_tend
#'
mean_trend_beta_vec <- function(form, beta, IP, timepoints){
    if(form %in% c("linear", "quadratic", "cubic", "L_down")){
        beta <- beta
    } else if(form == "M" || form == "W"){
        mu_IP1 <- as.numeric(t(beta) %*% c(0, IP[1]))
        beta_2 <- (0 - mu_IP1)/(IP[2] - IP[1])
        beta_3 <- (mu_IP1 - 0)/(IP[3] - IP[2])
        beta_4 <- (0 - mu_IP1)/(max(timepoints) - IP[3])
        beta <- c(beta, mu_IP1, beta_2, beta_3, beta_4)
    } else if(form == "L_up"){
        beta_slope <- beta
        beta_xintercept <- -beta_slope * IP
        beta <- c(beta_xintercept, beta_slope)
    }
    return(beta)
}

#' Inflection point check for \code{mean_trend}
#'
#' @param form character value specifying the type of time trend. Options
#'  include 'linear', 'quadratic', 'cubic', 'M', 'W', 'L_up', and 'L_down'.
#' @param beta vector specifying the appropriate parameters for functional
#' trend. See details of \code{\link{mean_trend}} for explanation for each form
#' @param IP vector specifying the inflection points. See details of
#' \code{\link{mean_trend}} for explanation for each form
#' @param timepoints numeric vector specifying the points to fit the functional
#'  trend.
#'
#' @return
#' Updated inflection point vector
IP_form_check <- function(form, beta, IP, timepoints){
    if (is.unsorted(IP)) {
        warning("Re-ordering IP from smallest to largest", call.=FALSE)
        IP <- IP[order(IP)]
    }
    if(form == "L_up"){
        if (length(IP) > 1) {
            stop("IP must be single value", call.=FALSE)
        }
        if (IP > max(timepoints) || IP < min(timepoints)) {
            stop("IP value of ", IP, " is outside of the range of timepoints\n
            [", round(min(timepoints), 2), ", ",round(max(timepoints), 2),"].
            Change either timepoints range or IP value.", call.=FALSE)
        }
    } else if (form == "M" || form == "W"){
        if (length(IP) != 3) {
            stop("IP must have three values", call.=FALSE)
        }
        if (any(IP == min(timepoints)) || any(IP == max(timepoints))) {
            stop("IP cannot equal min or max of timepoints", call.=FALSE)
        }
        if (sum(IP < max(timepoints)) != 3 || sum(IP > min(timepoints)) != 3) {
            warning("IP points outside of timepoint range", call.=FALSE)
        }
    } else if (form == "L_down"){
        IP <- -beta[1]/beta[2]
        if (IP > max(timepoints) || IP < min(timepoints)) {
            stop("IP value of ", IP, " is outside of the range of timepoints\n
                        [", round(min(timepoints), 2), ", ",
                round(max(timepoints),2), "].\n
                    IP is calculcated as -beta_yintercept/beta_slope.\n
                    Please try specifying these values again.",
                call.=FALSE)
        }
    }
    return(IP)
}
