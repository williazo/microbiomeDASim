.bioconduc_pkgs <- c("metagenomeSeq", "microbiomeDASim")
invisible(lapply(.bioconduc_pkgs, require, character.only=TRUE))
# loading in longitudinal microbiome data from Turnbaugh et. al 2009
data("mouseData")
mouseData

#aggregating the counts to the genus level
genus_mouseData <- metagenomeSeq::aggTax(mouseData, lvl="genus")
genus_mouseData

#restricting the
genus_mouseData <- filterData(genus_mouseData, present = 10, depth = 1000)
g_lnorm_mat <- MRcounts(genus_mouseData, norm=TRUE, log=TRUE)
set.seed(1234)
ex_feature <- sample(seq_len(dim(g_lnorm_mat)[1]), 1)
row.names(g_lnorm_mat)[ex_feature] #Anaerofustis
y <- g_lnorm_mat[ex_feature, ]

sim_feat <- data.frame(y, pData(genus_mouseData))
sim_feat <- sim_feat[order(sim_feat$mouseID, sim_feat$relativeTime), ]

sim_id_split <- split(sim_feat, sim_feat$mouseID)
sim_df <- lapply(sim_id_split, function(x){
    x$diet <- ifelse(any(x$diet=="Western"), "Western", x$diet)
    return(x)
})
sim_df <- data.frame(do.call(rbind, sim_df))
sim_df$diet <- as.factor(sim_df$diet)

with(sim_df, microbiomeDASim::ggplot_spaghetti(y, mouseID, relativeTime, group=diet))+
    xlab("Time")+
    ylab("Normalized Reads")+
    scale_color_discrete(name="Diet")+
    scale_linetype_discrete(name="Diet")+
    ggtitle(paste0("Genus: ", row.names(genus_mouseData[ex_feature, ])))+
    geom_vline(xintercept=21, col="black", lty=2)

ss_est <- fitTimeSeries(obj=genus_mouseData, formula=abundance~time*class,
                        feature=ex_feature, class="diet", id="mouseID",
                        time="relativeTime", norm=TRUE, log=TRUE, random=~1|id,
                        B=1000)
ss_est$timeIntervals
ggplot(data=ss_est$fit, aes(x=timePoints, y=fit))+
    geom_line()+
    geom_line(col="red", lty=2, aes(y=fit+1.96*se))+
    geom_line(col="red", lty=2, aes(y=fit-1.96*se))+
    geom_hline(col="black", yintercept=0)+
    geom_vline(xintercept=21, col="black", lty=2)

x_bar <- mean(sim_df$y)
s <- sd(sim_df$y)
gee_fit <- geepack::geeglm(y~1, id=mouseID, data=sim_df, corstr="ar1")
gee_fit <- geepack::geeglm(y~1, id=mouseID, data=sim_df, corstr="exch")

sim_mouse <- mvrnorm_sim_obs(id=sim_df$mouseID, time=sim_df$relativeTime,
                group=sim_df$diet, ref="BK", control_mean=x_bar, sigma=s,
                rho=0.03, corr_str="ar1", func_form="L_up", IP=21, beta=0.05)
with(sim_mouse$df, microbiomeDASim::ggplot_spaghetti(Y, ID, time, group=group))+
    xlab("Time")+
    ylab("Simulated Normalized Reads")+
    scale_color_discrete(name="Diet")+
    scale_linetype_discrete(name="Diet")+
    ggtitle(paste0("Genus: ", row.names(genus_mouseData[ex_feature, ])))+
    geom_vline(xintercept=21, col="black", lty=2)
