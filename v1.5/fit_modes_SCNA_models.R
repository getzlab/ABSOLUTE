## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.



fit_modes_SCNA_models = function( seg.obj, mode.tab, SCNA_model, mut.cn.dat, verbose=FALSE )
{
  Q = SCNA_model[["kQ"]]
  n.modes = nrow(mode.tab)
  obs <- seg.obj[["obs.scna"]]

  if (verbose) {
    print(paste("Optimizing SCNA_model | comb) for ",
                n.modes, " modes: ", sep=""))
  }
  
## TODO: absorb into SCNA_model
  log_ccf_dens = array(NA, dim=c(n.modes, length(obs[["W"]]), length(SCNA_model[["ccf_grid"]])) )
  dimnames(log_ccf_dens)[[3]] = SCNA_model[["ccf_grid"]]
  subclonal_scna_tab = array(NA, dim=c(n.modes, length(obs[["W"]]), 7) )
  dimnames(subclonal_scna_tab)[[3]] = c("CCF_hat", "CI95_low", "CI95_high", "subclonal_ix", "Pr_subclonal", "qs", "qc" )

  new_cols = c("genome mass", "sigma.h.hat", "theta.z.hat", "sigma.A.hat", "theta.Q.hat", "lambda.hat", "theta.0", "frac.het", "SCNA_LL", "entropy", "Kar_LL", "WGD", "combined_LL", "SSNV_LL", "SCNA_Theta_integral" )

  old_cols = colnames(mode.tab)
  mode.tab <- cbind(mode.tab, matrix(NA, nrow=nrow(mode.tab), ncol=length(new_cols))) 
  colnames(mode.tab) = c( old_cols, new_cols ) 
##
  mode_SCNA_models = list()

## optimize other parameters
  for (i in 1:n.modes) 
  {
    delta <- mode.tab[i, "delta"]
    b <- mode.tab[i, "b"]

    if(verbose) {
       print( paste("Optimizing PP mode #", i, sep=""))
    }

  ## optimize
    SCNA_model = SCNA_model_init(SCNA_model)
    SCNA_model = SCNA_model_calc_CCF_DP_loglik( obs, b, delta, SCNA_model, verbose=verbose )
    SCNA_model = calc_mode_seg_tabs( seg.obj, SCNA_model, b, delta )

## Annotate SCNA clonality summary for SSNV models
    res = get_subclonal_SCNA_info( obs, b, delta, SCNA_model, mut.cn.dat )
    subclonal_scna_tab[i,,] = res[["subclonal_scna_tab"]]
    log_ccf_dens[i,,] = res[["log_ccf_dens"]]
##
    mode.tab = fill_mode.tab_row( seg.obj, SCNA_model, obs, b, delta, mode.tab, i )

    mode_SCNA_models[[i]] = SCNA_model
  }

  subclonal_SCNA_res = list(subclonal_SCNA_tab = subclonal_scna_tab, log_CCF_dens = log_ccf_dens)

  return( list(mode.tab=mode.tab, mode_SCNA_models=mode_SCNA_models, subclonal_SCNA_res=subclonal_SCNA_res, mode.flag=NA ) )
}




WeighSampleModes <- function(mode.res) 
{
## combined various scores
  mode.tab = mode.res[["mode.tab"]]

## Only using SCNA score for now..
  mode.res[["mode.tab"]][, "combined_LL"] = mode.tab[,"SCNA_LL"] # + mode.tab[,"Kar_LL"] + mode.tab[,"SSNV_LL"] 

  LL = mode.res[["mode.tab"]][, "combined_LL"]
  if( !all( is.finite(LL)) ) { stop("Non-finite mode combined_LL!") }

  dens = exp(LL - LogAdd(LL))
  mode.res[["mode.tab"]] <- cbind(mode.res[["mode.tab"]], dens)
  ix = order(LL, decreasing=TRUE )
  
  mode.res <- ReorderModeRes(mode.res, ix)
  
  return(mode.res)
}

fill_mode.tab_row = function( seg.obj, SCNA_model, obs, b, delta, mode.tab, i )
{
   Q = SCNA_model[["kQ"]]
## This is leftover from ancient code - copies fields from fit SCNA_model into mode.tab.   Would it be easier to just save a list of SCNA_model objects??
    Theta_hat = SCNA_model_Theta_tx(SCNA_model[["Theta"]])  ## recover natural units
    mode.tab[i, "theta.0"] <- SCNA_model[["theta.0"]]
    mode.tab[i, "theta.Q.hat"] <- Theta_hat["theta.Q"]
    mode.tab[i, "sigma.A.hat"] <- Theta_hat["sigma.A"]
    mode.tab[i, "sigma.h.hat"] <- SCNA_model[["sigma.h"]]
    mode.tab[i, "lambda.hat"] <- Theta_hat["lambda"]
    mode.tab[i, "SCNA_LL"] <- SCNA_model[["LL"]]
#    mode.tab[i, "theta.z.hat"] <-  sum(SCNA_model[["mix.w"]][c("unif","exp")])
    mode.tab[i, "theta.z.hat"] = sum( SCNA_model[["seg.qz.tab"]][,(Q + 1)] * obs$W)
    ## compute % non-clonal genome
    frac.het = sum(obs[["W"]] * SCNA_model[["seg.z.tab"]])  
    if( !is.finite(frac.het)) { stop() }

    mode.tab[i, "frac.het"] <- frac.het
    ## calculate allelic-balance
    if( seg.obj[["copy_num_type"]] == "allelic" )  {
      mode.tab[i,"genome mass"] <- 2 * sum(c((1:Q)-1) * colSums( SCNA_model[["seg.q.tab"]] * obs[["W"]]))
    }
    if( seg.obj[["copy_num_type"]] == "total" ) {
      mode.tab[i,"genome mass"] <- 1 * sum(c((1:Q)-1) * colSums( SCNA_model[["seg.q.tab"]] * obs[["W"]]))
    }

    mode.tab[i, "WGD"] = SCNA_model[["WGD"]]

    ## weighted entropy average over segs
    mode.tab[i, "entropy"] <- CalcFitEntropy(obs, SCNA_model[["seg.qz.tab"]] )

    return(mode.tab)
}

ReorderModeRes <- function(mode.res, ix, DROP=FALSE) 
{
   mode.res[["mode_SCNA_models"]] = mode.res[["mode_SCNA_models"]][ix] 
   mode.res[["mode.tab"]] <- mode.res[["mode.tab"]][ix,, drop=DROP]
   mode.res[["mode.clust.p"]] <- mode.res[["mode.clust.p"]][ix , , drop=DROP]  ## From Kar model

# ? does this crash with no MAF??
   mode.res[["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]] = mode.res[["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]][ix, , , drop=DROP]
   mode.res[["subclonal_SCNA_res"]][["log_CCF_dens"]] = mode.res[["subclonal_SCNA_res"]][["log_CCF_dens"]][ix, , , drop=DROP]

   ## only exists if MAF supplied.
   if (!is.null(mode.res[["modeled.muts"]])) 
   {
      mode.res[["SSNV.ccf.dens"]] = mode.res[["SSNV.ccf.dens"]][ix,,,drop=DROP]
      mode.res[["modeled.muts"]] <- mode.res[["modeled.muts"]][ix]
      mode.res[["mode_SSNV_models"]] = mode.res[["mode_SSNV_models"]][ix]
   }
   mode.res[["mode.posts"]] <- mode.res[["mode.posts"]][ix]

   return(mode.res)
}


## dont call non-aneuploid if mut_dat is present
GetCallStatus <- function(mode.res, seg.w) {
  status <- "called"
  q.tab <- mode.res[["mode_SCNA_models"]][[1]][["seg.qz.tab"]]
  q.tab <- q.tab[, c(1:  (ncol(q.tab)) )]
  
  Q <- ncol(q.tab) - 1
  max.q <- (apply(q.tab, 1, which.max))

  peak_masses <- rep(0, Q)
  for (i in 1:Q) {
    ix <- which(max.q==i) 
    peak_masses[i] <- sum(seg.w[ix])
  }
  
  b <- mode.res[["mode.tab"]][1,"b"]
  
  ## don't count 0-CN state if b is too small - could be due to IBD, not LOH
  if (peak_masses[1] < 0.01 & b < 0.15) {
    peak_masses[1] <- 0
  }

  six <- order(peak_masses, decreasing=TRUE)

  if (peak_masses[six[3]] < 0.0001) {
    if (peak_masses[six[2]] < 0.0001) { 
      if (mode.res[["mode.tab"]][1,"sigma.h.hat"] < 0.02 &
          is.null(mode.res[["modeled.muts"]])) { 
        status <- "non-aneuploid" 
      }
      if (mode.res[["mode.tab"]][1,"sigma.h.hat"] >= 0.02 ) {
        status <- "low purity"
      }
    }
  }

   if (mode.res[["mode.tab"]][1, "entropy"] > 0.2) {
     status <- "high entropy"
   }

   if (mode.res[["mode.tab"]][1, "frac.het"] > 0.2) {
     status <- "high non-clonal"
   }

   return(status)
}
