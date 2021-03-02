## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


## The purpose of this function is to add additional purity/ploidy modes to the 
## initial set of provisional modes.   This optimization assumes genome_mass=2 
## and searches for the optimal alpha by fitting the somatic mutations.

old_run_diploid_snv_purity_opt = function(obs_scna, mut.cn.dat, SSNV_model,
				 dom, d_res = 0.1, verbose=FALSE) 
{
  seg = obs_scna$segtab
  ## expected copy-number, should be 1.0
  e_cr = sum(seg[, "W"] * seg[, "copy_num"])  
  
  ## remove muts that are not on diploid regs
  scna.ix = which(abs(seg[, "copy_num"] - e_cr) > 0.1)
  
  if (length(scna.ix) > 0) {
    ##  works for HSCR or total CR
    if ("A1.ix" %in% colnames(mut.cn.dat)) {
      nix = mut.cn.dat[, "A1.ix"] %in% scna.ix | mut.cn.dat[,"A2.ix"] %in% scna.ix
    } else {
      nix = mut.cn.dat[, "mut_seg_ix"] %in% scna.ix 
    }     
    
    mut.cn.dat = mut.cn.dat[!nix, , drop=FALSE]
  }
  
  if (nrow(mut.cn.dat) == 0) {
    if (verbose) {
      print("No muts left on unaltered genome")
    }
    return(NA)
  }
  
  comb_1d_ll = function(par, mut.cn.dat, dom) {
    alpha = par
    
    ## FIXME: This looks like it should be a ||
    if (alpha < dom[1] | alpha > dom[2]) {
      return( .Machine$double.xmax) 
    } 
    
    mode_info = alpha
    names(mode_info)="alpha"

    mult_res = fit_SSNV_model(mut.cn.dat, mode_info, SSNV_model, print_fit=FALSE)
    
    ll = sum(mult_res$post_Prs[, "LL"], na.rm=TRUE)  + dDirichlet(mult_res$som_theta_Q_MAP, 
                                                                  SSNV_model[["kPiSomThetaQ"]], log.p=TRUE)
    return(-ll)
  }
  
  mut.cn.dat = cbind(mut.cn.dat, q_hat=2, HS_q_hat_1=1, HS_q_hat_2=1)
  
  d_grid = seq(dom[1], dom[2], length=100)
  ll_grid = rep(NA, length(d_grid))
  for (i in seq_along(d_grid)) {
    ll_grid[i] = -comb_1d_ll(d_grid[i], mut.cn.dat, dom)
    cat("~")
  }
  cat("\n")

  d_ll = diff(ll_grid) 
  sd_ll = d_ll > 0
  rs = rbind(sd_ll, c(0, sd_ll[1:(length(sd_ll - 1))]))
  zc = which(rs[1, ]== 0 & rs[2, ]==1)
  mode_ix = zc
  if (ll_grid[length(ll_grid)] > ll_grid[length(ll_grid) - 1]) { 
    mode_ix = c(mode_ix, length(ll_grid)) 
  }
  if (ll_grid[1] > ll_grid[2]) { 
    mode_ix = c(mode_ix, 1) 
  }
  mode_tab = array( NA, dim=c(length(mode_ix), 3))
  
  a = d_grid[mode_ix]
  tau = 2
  # solve so comb[1] == e_cr
  tau_p = (2 * (1 - a) + 2 * a) / (e_cr * a)  -  2 * (1 - a) / a
  
  ##  convert to tau corresponding to ploidy=2 for this sample
  b = 2 * (1 - a) / (2 * (1 - a) + a * tau_p)
  delta = a / (2 * (1 - a) + a * tau_p)
  
  mode_tab[, 1] = b
  mode_tab[, 2] = log(delta)
  mode_tab[, 3] = ll_grid[mode_ix]
     
  return(mode_tab)
}






run_diploid_snv_purity_opt = function(obs_scna, mut.cn.dat, SSNV_model,
				 dom, d_res = 0.1, verbose=FALSE) 
{
  ## expected copy-number, should be 1.0
  e_cr = obs_scna[["e.cr"]]
  
  ## remove muts that are not on diploid regs
  scna.ix = which(abs(obs_scna[["d"]] - e_cr) > 0.1)
  
  if (length(scna.ix) > 0) {
    ##  works for HSCR or total CR
    if ("A1.ix" %in% colnames(mut.cn.dat)) {
      nix = mut.cn.dat[, "A1.ix"] %in% scna.ix | mut.cn.dat[,"A2.ix"] %in% scna.ix
    } else {
      nix = mut.cn.dat[, "mut_seg_ix"] %in% scna.ix 
    }     
    
    mut.cn.dat = mut.cn.dat[!nix, , drop=FALSE]
    SSNV_model[["d_f_grid"]] = SSNV_model[["d_f_grid"]][!nix,] # this mod does not survive
  }
  
  if (nrow(mut.cn.dat) == 0) {
    if (verbose) {
      print("No muts left on unaltered genome")
    }
    return(NA)
  }
  
  comb_1d_ll = function(par, mut.cn.dat, dom) {
    alpha = par
    
    if (alpha < dom[1] ) { alpha=dom[1] }
    if( alpha > dom[2] ) { alpha=dom[2] }
    
    mode_info = alpha
    names(mode_info)="alpha"

    mult_res = fit_SSNV_model(mut.cn.dat, mode_info, SSNV_model, print_fit=FALSE)
    
    ll = sum(mult_res$post_Prs[, "LL"], na.rm=TRUE)  + dDirichlet(mult_res$som_theta_Q_MAP, 
                                                                  SSNV_model[["kPiSomThetaQ"]] , log.p=TRUE)
    return(-ll)
  }
  
  mut.cn.dat = cbind(mut.cn.dat, q_hat=2, HS_q_hat_1=1, HS_q_hat_2=1)
  
  d_grid = seq(dom[1], dom[2], length=25)
  opt_alpha = array(NA, dim=c( length(d_grid), 2))

  for (i in seq_along(d_grid))
  {
    opt = nlm(f=comb_1d_ll, p=d_grid[i], mut.cn.dat=mut.cn.dat, dom=dom)
#    mode_tab[i, ] = c(0, opt$estimate, -opt$minimum)
    opt_alpha[i, ] = c(opt$estimate, -opt$minimum)

    cat("~")
  }
  cat("\n")


  nix = opt_alpha[,1] < dom[1] | opt_alpha[,1] > dom[2] 
  opt_alpha = opt_alpha[!nix,]

  a = opt_alpha[,1]
  tau = 2
  # solve so comb[1] == e_cr
  tau_p = (2 * (1 - a) + 2 * a) / (e_cr * a)  -  2 * (1 - a) / a
  
  ##  convert to tau corresponding to ploidy=2 for this sample
  b = 2 * (1 - a) / (2 * (1 - a) + a * tau_p)
  delta = a / (2 * (1 - a) + a * tau_p)
  
  mode_tab = array(NA, dim=c( length(b), 3))

  mode_tab[, 1] = b
  mode_tab[, 2] = log(delta)
  mode_tab[, 3] = opt_alpha[,2]
     
  return(mode_tab)
}




