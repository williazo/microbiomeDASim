#' Convert simulated output to MRexperiment object
#'
#' In order to allow investigators to more easily incorporate simulated data,
#' this package converts the raw output into an MRexperiment object used in the
#' \code{\link{metagenomeSeq}} package.
#'
#' @param obj output from either \code{\link{gen_norm_microbiome}} or
#' \code{\link{mvrnorm_sim}}
#' @param missing logical indicator for objects from \code{\link{mvrnorm_sim}}.
#' If missing = TRUE then create MRexperiment object with \code{Y_obs} else
#' use \code{Y}.
#'
#' @examples
#' bug_gen <- gen_norm_microbiome(features=6, diff_abun_features=3,
#'                                n_control=30, n_treat=20, control_mean=2,
#'                                sigma=2, num_timepoints=4, t_interval=c(0, 3),
#'                                rho=0.9, corr_str="compound", func_form="M",
#'                                beta=c(4, 3), IP=c(2, 3.3, 6),
#'                                missing_pct=0.2, missing_per_subject=2,
#'                                miss_val=0, asynch_time=TRUE)
#' bug_gen_MR <- simulate2MRexperiment(bug_gen)
#' class(bug_gen_MR)
#'
#' @export
#'
#' @return
#' An MRexperiment object
#'
#' @importFrom metagenomeSeq newMRexperiment
#' @importFrom Biobase AnnotatedDataFrame
simulate2MRexperiment <- function(obj, missing=FALSE){
    #obj comes from gen_norm_microbiome function
    if(length(obj)==2){
        obj_pheno <- Biobase::AnnotatedDataFrame(obj$bug_feat)
        obj <- metagenomeSeq::newMRexperiment(counts=obj$Y, phenoData=obj_pheno)
    }
    #obj from mvrnorm_sim
    else if(length(obj)==7 && missing==FALSE){
    obj_pheno <- Biobase::AnnotatedDataFrame(obj$df[, c("ID","time", "group")])
    obj <- metagenomeSeq::newMRexperiment(counts=t(obj$Y), phenoData=obj_pheno)
    }
    else if(length(obj)==7 && missing==TRUE){
    obj_pheno <- Biobase::AnnotatedDataFrame(obj$df[, c("ID","time", "group")])
    obj <- metagenomeSeq::newMRexperiment(counts=t(obj$Y_obs),
                                            phenoData=obj_pheno)
    }
    else{
        stop("unrecognized object. See details for accepted objects",
                call.=FALSE)
    }
    return(obj)
}

#' Convert simulated output to phyloseq object
#'
#' This function will convert simulated data into a \code{\link{phyloseq}}
#' object.
#'
#' @param obj output from either \code{\link{gen_norm_microbiome}} or
#' \code{\link{mvrnorm_sim}}
#' @param missing logical indicator for objects from \code{\link{mvrnorm_sim}}.
#' If missing = TRUE then create MRexperiment object with \code{Y_obs} else
#' use \code{Y}.
#'
#' @examples
#' bug_gen <- gen_norm_microbiome(features=6, diff_abun_features=3,
#'                                n_control=30, n_treat=20, control_mean=2,
#'                                sigma=2, num_timepoints=4, t_interval=c(0, 3),
#'                                rho=0.9, corr_str="compound", func_form="M",
#'                                beta=c(4, 3), IP=c(2, 3.3, 6),
#'                                missing_pct=0.2, missing_per_subject=2,
#'                                miss_val=0, asynch_time=TRUE)
#' bug_gen_phyloseq <- simulate2MRexperiment(bug_gen)
#' class(bug_gen_phyloseq)
#'
#' @importFrom phyloseq phyloseq otu_table sample_data
#'
#' @export
#'
#' @return
#' A phyloseq object
simulate2phyloseq <- function(obj, missing=FALSE){
    #obj comes from gen_norm_microbiome function
    if(length(obj)==2){
        obj_tbl <- phyloseq::otu_table(obj$Y, taxa_are_rows=TRUE)
        obj_pheno <- phyloseq::sample_data(obj$bug_feat)
        obj <- phyloseq::phyloseq(obj_tbl, obj_pheno)
    }
    #obj from mvrnorm_sim
    else if(length(obj)==7 && missing==FALSE){
        obj_tbl <- phyloseq::otu_table(t(obj$Y), taxa_are_rows=TRUE)
        obj_pheno <- phyloseq::sample_data(obj$df[, c("ID","time", "group")])
        obj <- phyloseq::phyloseq(obj_tbl, obj_pheno)
    }
    else if(length(obj)==7 && missing==TRUE){
        obj_tbl <- phyloseq::otu_table(t(obj$Y_obs), taxa_are_rows=TRUE)
        obj_pheno <- phyloseq::sample_data(obj$df[, c("ID","time", "group")])
        obj <- phyloseq::phyloseq(obj_tbl, obj_pheno)
    }
    else{
        stop("unrecognized object. See details for accepted objects",
                call.=FALSE)
    }
    return(obj)
}
