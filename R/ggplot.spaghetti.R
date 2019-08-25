#'Spaghetti Plots using \code{ggplot2}
#'
#'This function allows the user to create spaghetti plots for individuals
#' with time varying covariates. You can also break this down into subgroups
#' to analyze different trentds.
#'
#'Note that the data must be in long format.
#'
#'@param y This is the y-axis parameter to specify.
#'Generally it is a continuous variable.
#'@param id This is the id parameter that identifies the unique individuals or
#' units.
#'@param time This is the time vector and must be numeric.
#'@param alpha Scalar value between \[0,1\] that specifies the transparencey
#' of the lineplots.
#'@param method Character value that specifies which type of method to use for
#' fitting. Optional methods come from \code{\link[ggplot2]{stat_smooth}}
#' function.
#'@param jit Scalar value that specifies how much you want to jitter each
#' individual observation. Useful if many of the values share the same y values
#' at a time point.
#'@param group Specifies a grouping variable to be used, and will plot it by
#' color on one single plot.
#'
#'@return Plots a time series data by each individual/unit with group trends
#' overlayed.
#'
#'@examples
#'library(ggplot2)
#'num_subjects_per_group <- 15
#'sim_obj <- mvrnorm_sim(n_control=num_subjects_per_group,
#'                        n_treat=num_subjects_per_group,
#'                        control_mean=5, sigma=1, num_timepoints=5,
#'                        rho=0.95, corr_str='ar1', func_form='linear',
#'                        beta=c(0, 0.25),
#'                        missing_pct=0.6, missing_per_subject=2)
#'
#'with(sim_obj$df, suppressWarnings(ggplot_spaghetti(y=Y_obs, id=ID, time=time,
#'                                                   jit=0.1, group=group)))+
#'   labs(title="Simulated Microbiome Data from Multivariate Normal",
#'        y="Normalized Reads", x="Time") +
#'   scale_linetype_manual(values=c("solid","dashed"), name="Group") +
#'   scale_color_manual(values=c("#F8766D", "#00BFC4"), name="Group")
#'
#'@import ggplot2
#'@importFrom stats runif
#'
#'@export
#'
ggplot_spaghetti <- function(y, id, time, alpha = 0.2, method = "loess",
                                jit = 0.0, group = NULL){
    fact <- ifelse(jit == 0, 0, 1)
    if(is.factor(group) == FALSE) group <- factor(group)
    gg_dat <- data.frame(y, id, time, group)
    gg_dat <- subset(gg_dat, !is.na(y))
    ids <- as.character(unique(gg_dat$id))
    groups <- unique(as.character(gg_dat$group))
    if(length(groups)>13){
        lty_group <- as.factor(rep("solid", nrow(gg_dat)))
    } else {
        lty_group <- gg_dat$group
    }
    gg_dat$lty_group <- lty_group
    base <- ggplot() + xlab("") + ylab("")
    for(i in ids){
        for (j in groups){
            ry <- runif(1, min = 0, max = jit)
            rx <- runif(1, min = 0, max = jit)
            gg_dat_grp <- gg_dat[gg_dat$id == i & gg_dat$group == j,]
            gg_dat_grp$time <- jitter(gg_dat_grp$time, factor = fact,
                                        amount = rx)
            gg_dat_grp$y <- jitter(gg_dat_grp$y, factor = fact, amount = ry)
            if(nrow(gg_dat_grp) >= 1){
                base <- base +
                        geom_point(data = gg_dat_grp,
                                aes(x = time, y = y, col = group,
                                    linetype = lty_group), alpha = alpha)
            }
            if(nrow(gg_dat_grp) >= 2){
                base <- base +
                    geom_line(data = gg_dat_grp,
                                aes(x = time, y = y, col = group,
                                    linetype = lty_group), alpha = alpha)
            }
        }
    }
    base <- base +
        stat_smooth(data = gg_dat,
                    aes(x = time, y = y, group = group, col = group,
                            linetype = lty_group),
                    lwd = 2.5, method = method, se = FALSE)
    return(base)
}

