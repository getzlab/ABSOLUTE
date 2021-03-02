## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

calc_muts_ccf_postrior_LL_grid = function( mut.cn.dat, alpha, SSNV_model )
{ 
## TODO: add X-chr support
  ccf_grid = SSNV_model[["ccf_grid"]]

## Calc posterior density of CCF | d
  ccf.LL = matrix(NA, nrow=nrow(mut.cn.dat), ncol=length(ccf_grid))
  for (i in seq_len(nrow(mut.cn.dat))) 
  {
    if (mut.cn.dat[i, "q_hat"] == 0) { 
      next 
    }
#    ccf.dens[i, ] = calc_ccf_posterior_dens_grid(mut.cn.dat[i, "alt"], mut.cn.dat[i, "ref"], SSNV_model, alpha, mut.cn.dat[i, "q_hat"] )

    ccf.LL[i, ] = calc_ccf_posterior_LL_grid(mut.cn.dat[i, "alt"], mut.cn.dat[i, "ref"], SSNV_model, alpha, mut.cn.dat[i, "q_hat"] )
  }
##
  return(ccf.LL)
}


fit_SSNV_model = function(mut.cn.dat, mode_info, SSNV_model, allelic_subclonal_scna_tab, scna_log_ccf_dens, tol=0.001, max.iter=25, verbose=FALSE )
{
  clonal_scna_mut_ix = mut.cn.dat[,"clonal_scna_mut_ix"]

 ## assume clonal SCNA initially
  post.ccf.LL.grid = calc_muts_ccf_postrior_LL_grid( mut.cn.dat, mode_info["alpha"], SSNV_model )
  ssnv.ccf.dens = exp( post.ccf.LL.grid - LogAdd( post.ccf.LL.grid))

  iter=1
  cur.loglik = -Inf

  while( 1 )
  {
    post_Prs = NA
# SSNVs on clonal SNCAs
    if( sum( clonal_scna_mut_ix ) > 0 )
    {
      clonal_SCNA_res = eval_SNV_models_evidence(mut.cn.dat[ clonal_scna_mut_ix,, drop=FALSE ], mode_info, post.ccf.LL.grid[ clonal_scna_mut_ix,, drop=FALSE], SSNV_model )

      mut.ev.mat = matrix( NA, nrow=nrow(mut.cn.dat), ncol=ncol(clonal_SCNA_res[["mut.ev.mat"]]) )
      colnames(mut.ev.mat) = colnames(clonal_SCNA_res[["mut.ev.mat"]])

      post_Prs = matrix( NA, nrow=nrow(mut.cn.dat), ncol=ncol(clonal_SCNA_res[["post_Prs"]]))
      colnames(post_Prs) = colnames(clonal_SCNA_res[["post_Prs"]])
      post_Prs[clonal_scna_mut_ix,] = clonal_SCNA_res[["post_Prs"]]
    }
    else
    {
      clonal_SCNA_res = list( "mut.ev.mat"=NA, "som_mut_Q_tab"=NA)
    }

# SSNVs on subclonal SCNAs
    if( sum( !clonal_scna_mut_ix ) > 0 )
    {
      SSNV_on_subclonal_SCNA_res = allelic_calc_sample_muts_on_subclonal_scna( mut.cn.dat[ !clonal_scna_mut_ix,, drop=FALSE ], mode_info, allelic_subclonal_scna_tab, scna_log_ccf_dens, SSNV_model ) 

      mut.ev.mat = matrix( NA, nrow=nrow(mut.cn.dat), ncol=ncol(SSNV_on_subclonal_SCNA_res[["mut.ev.mat"]]) )
      colnames(mut.ev.mat) = colnames(SSNV_on_subclonal_SCNA_res[["mut.ev.mat"]])
      if(sum(clonal_scna_mut_ix) == 0 ) ## in case no muts on clonal SCNA..
      {
        post_Prs = matrix( NA, nrow=nrow(mut.cn.dat), ncol=ncol(SSNV_on_subclonal_SCNA_res[["post_Prs"]]) )
        colnames(post_Prs) = colnames(SSNV_on_subclonal_SCNA_res[["post_Prs"]])
      }
      post_Prs[!clonal_scna_mut_ix,] = SSNV_on_subclonal_SCNA_res[["post_Prs"]]
    }  
    else
    {
      H.ev = matrix( NA, nrow=0, ncol=4 ) 
      colnames(H.ev) = c("H1", "H2", "H3", "H4")
      SSNV_on_subclonal_SCNA_res = list( "mut.ev.mat"=NA, "ssnv.ccf.dens"=NA, "H.ev"=H.ev, "som_mut_Q_tab"=NA)
    }
    SSNV_model[["SSNV_on_subclonal_SCNA_res"]] = SSNV_on_subclonal_SCNA_res

## put them together
    som_mut_Q_tab = matrix( NA, nrow=nrow(mut.cn.dat), ncol=SSNV_model[["kQ"]] )
    som_mut_Q_tab[clonal_scna_mut_ix,] = clonal_SCNA_res[["som_mut_Q_tab"]] 
    som_mut_Q_tab[!clonal_scna_mut_ix,] = SSNV_on_subclonal_SCNA_res[["som_mut_Q_tab"]] 

    loglik = sum(post_Prs[,"LL"])
    cond <- abs(cur.loglik - loglik) / abs(cur.loglik)

    print(paste("loglik = ", round(loglik,4), sep=""))
    print(paste("cond = ", round(cond, 4), sep=""))
 
    if (verbose) { 
      print(SSNV_model[["mut_class_w"]]) 
    }
    if (verbose) {  
      cat("som_theta_Q_MAP: ")
      print(round( SSNV_model[["som_theta_Q_mode"]], 5)) 
    }

    if(( iter > 1 & cond < tol) || (iter >= max.iter)) { break }
#    break  # turn off Hierarchical fitting in a single sample

    iter = iter+1
    cur.loglik = loglik

    SSNV_model = update_SSNV_mixture_weights( SSNV_model, som_mut_Q_tab, post_Prs ) 
  }   

## update ccf dens for ssnvs on subclonal scnas
  ssnv.ccf.dens[!clonal_scna_mut_ix,] = SSNV_on_subclonal_SCNA_res[["ssnv.ccf.dens"]]

  mut.ev.mat[clonal_scna_mut_ix,] = clonal_SCNA_res[["mut.ev.mat"]]
  mut.ev.mat[!clonal_scna_mut_ix,] = SSNV_on_subclonal_SCNA_res[["mut.ev.mat"]]

## only exists for SSNVs on subclonal SCNAs
  H.ev = matrix( NA, nrow=nrow(mut.cn.dat), ncol=ncol(SSNV_on_subclonal_SCNA_res[["H.ev"]]) )
  colnames(H.ev) = colnames(SSNV_on_subclonal_SCNA_res[["H.ev"]])
  H.ev[!clonal_scna_mut_ix,] = SSNV_on_subclonal_SCNA_res[["H.ev"]]

  if( any( is.nan(ssnv.ccf.dens))) { stop("NaN in ssnv.ccf.dens") }

  if( nrow(ssnv.ccf.dens) != nrow(post_Prs) ) { stop() }

  return( list( "post_Prs" = post_Prs, "som_theta_Q_MAP" = SSNV_model[["som_theta_Q_mode"]], "ssnv.ccf.dens" = ssnv.ccf.dens, "som_mut_Q_tab" = som_mut_Q_tab, "mut.ev.mat" = mut.ev.mat, "H.ev"=H.ev, "SSNV_model"=SSNV_model ) )
}

## Not clear this is optimal - better to summarize    SSNV_model directly..
## Computes CCF 95 CI for output to MAF tab
calc_CCF_95CI = function(mut_dat, ccf_dens, mode_info, modal_q_s, pr_somatic_clonal, SSNV_model) 
{
  ccf_grid = SSNV_model[["ccf_grid"]] 
  muts = mut_dat
  alpha = mode_info["alpha"]
  cov = rowSums(muts[, c("alt", "ref")])
  fhat = muts[, "alt"] / cov
  Q = muts[, "q_hat"]
  
  q_s = rep(NA, length(pr_somatic_clonal))
  c.ix = pr_somatic_clonal > 0.5
  q_s[c.ix] = modal_q_s[c.ix]
  q_s[!c.ix] = 1 
  
  som_delta = alpha / (2 * (1 - alpha) + alpha * Q)
#  f_skew =  SSNV_model[["SSNV_skew"]]

  f_s = q_s * som_delta    
  mut_scale = 1 / f_s   
  cell_frac = fhat * mut_scale  # * f_skew^-1
#  cell_mult = fhat / som_delta  # * f_skew^-1
  
  ccf_ci95 = matrix(NA, nrow=nrow(ccf_dens), ncol=2)
  ccf_hat = rep(NA, nrow(mut_dat))
  
  for (i in seq_len(nrow(mut_dat))) 
  {
    if (mut_dat[i, "q_hat"] == 0 | any(is.nan(ccf_dens[i,])) ) { 
      next 
    }
    
    ccf_hat[i] = ccf_grid[which.max(ccf_dens[i, ])]
    
    ecdf = cumsum(ccf_dens[i, ])
    ccf_ci95[i, ] = approx(x=ecdf, y=ccf_grid, xout=c(0.025, 0.975))$y
  }
  nix1 = is.na(ccf_ci95[, 1])
  ccf_ci95[nix1, 1] = min(ccf_grid)
  nix2 = is.na(ccf_ci95[,2])
  ccf_ci95[nix2, 2] = max(ccf_grid)
  
  ## Round up in last bin.   TODO round down in 1st bin 
  ix = ccf_ci95[, 2] > ccf_grid[length(ccf_grid) - 1]
  ccf_ci95[ix, 2] = 1.0
  
  ix = mut_dat[i, "q_hat"] == 0
  ccf_ci95[ix, ] = NA
  
  res = cbind( ccf_hat, ccf_ci95)
  colnames(res) = c("ccf_hat", "ccf_CI95_low", "ccf_CI95_high")
  
  return(res)
}



## These fields are output in the ABS MAF for each SSNV..
## Does not take subclonal SCNA into account
## ??  Is this still useful ?
ClassifySomaticVariants <- function(prs, pr.thresh)
{
  subclonal.ix <- prs[, "Pr_subclonal"] > pr.thresh
  subclonal.wt0.ix <- prs[, "Pr_subclonal_wt0"] > pr.thresh
  clonal.ix <- prs[, "Pr_somatic_clonal"] > pr.thresh
  wt0.ix <- prs[, "Pr_wt0"] > pr.thresh
  ge2.ix <- prs[, "Pr_ge2"] > pr.thresh

  clonal.het.ix <- prs[, "Pr_somatic_clonal"] * (1 - prs[, "Pr_wt0"]) > pr.thresh 
  homozygous.ix <- rowSums(prs[, c("Pr_wt0", "Pr_subclonal_wt0"), drop=FALSE]) > pr.thresh
    
  res <- cbind(subclonal.ix, subclonal.wt0.ix, clonal.ix, wt0.ix, clonal.het.ix, 
               ge2.ix, homozygous.ix)
  colnames(res) <- c("subclonal.ix", "subclonal_wt0.ix", "clonal.ix", "wt0.ix", 
                     "clonal_het.ix", "ge2.ix", "homozygous.ix")
  
  return(res)
}


get_SSNV_on_clonal_CN_multiplicity_densities = function( seg.dat, mut.dat, af_post_pr, grid_mat, verbose=FALSE )
{
  ## remove muts on HZdels 
  hz.del.flag <- mut.dat[, "q_hat"] == 0
  nix = hz.del.flag

  alpha <- mut.dat[1, "purity"]
  Q <- mut.dat[, "q_hat"]
  som.delta <- alpha / (2 * (1 - alpha) + alpha * Q)
#  som.delta = som.delta * SSNV_skew^-1

  mult_grid = matrix( NA, nrow=nrow(grid_mat), ncol=ncol(grid_mat))
  mult_dens = matrix( NA, nrow=nrow(af_post_pr), ncol=ncol(af_post_pr))

  mult_grid[!nix, ] = (grid_mat / som.delta)   [!nix,, drop=FALSE]
  mult_dens[!nix, ] = (af_post_pr * som.delta) [ !nix,, drop=FALSE]

#  bad_rows = which( apply( is.nan(mult_dens), 1, sum ) > 0 )
#  if( any(is.nan(mult_dens)) ) { stop() }

  return( list("mult_dens"=mult_dens, "mult_grid"=mult_grid))
}


get_SSNV_on_subclonal_CN_multiplicity_densities = function( seg.dat, mut.dat, af_post_pr, grid_mat, verbose=verbose )
{
   stop()
}

