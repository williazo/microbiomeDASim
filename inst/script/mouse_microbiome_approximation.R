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
ex_feature <- sample(seq_len(dim(g_lnorm_mat)[1]), 1)
y <- g_lnorm_mat[ex_feature, ]

test_df <- data.frame(y, pData(genus_mouseData))
test_df <- test_df[order(test_df$mouseID, test_df$relativeTime), ]

test <- split(test_df, test_df$mouseID)
test <- lapply(test, function(x){
    x$diet <- ifelse(any(x$diet=="Western"), "Western", x$diet)
    return(x)
})
test <- data.frame(do.call(rbind, test))

with(test, microbiomeDASim::ggplot_spaghetti(y, mouseID, relativeTime, group=diet))+
    xlab("Time")+
    ylab("")+
    scale_color_discrete(name="Diet")+
    scale_linetype_discrete(name="Diet")+
    ggtitle(paste0("Genus: ", row.names(genus_mouseData[ex_feature, ])))+
    geom_vline(xintercept=21, col="black", lty=2)

gen_norm_microbiome()
