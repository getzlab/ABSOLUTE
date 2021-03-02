## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


init_SSNV_model = function( kQ, SSNV_skew, N_mut )
{ 
   SSNV_model = list()
 # Prior on somatic mutation multiplicity values [1...Q]
   SSNV_model[["kPiSomThetaQ"]] = c(250, 50, rep(2, (kQ - 2)))
   SSNV_model[["kPiSomThetaQ"]] = SSNV_model[["kPiSomThetaQ"]] / sum(SSNV_model[["kPiSomThetaQ"]]) * N_mut * 3  # prior is 3x stronger than data  
   SSNV_model[["kPiSomThetaQ"]][  SSNV_model[["kPiSomThetaQ"]] < 2 ] = 2

   SSNV_model[["mut_class_w"]] = list(SM = 0.5, GL = 0, SC = 0.5, OL = 1e-3, Pi_SM = 15, Pi_SC = 15)
   SSNV_model[["kQ"]] = kQ

   SSNV_model[["d_ccf"]] = 0.01
   SSNV_model[["ccf_grid"]] = seq(0, 1, by=0.01)

   SSNV_model[["SSNV.on.CN0"]] = 1e-10

   SSNV_model[["lambda"]] = 25
   SSNV_model[["W_U"]] = 0.25 
   SSNV_model[["W_E"]] = 0.75

### Initialize dir modes
  pi_som_theta_q = SSNV_model[["kPiSomThetaQ"]]
  mut_class_w = SSNV_model[["mut_class_w"]]

  ## Dir mode
  SSNV_model[["som_theta_Q_mode"]] = (pi_som_theta_q - 1) / (sum(pi_som_theta_q) - length(pi_som_theta_q)) 
  
  ## mut_class_w Dir mode
  mcv = c(mut_class_w$Pi_SM, mut_class_w$Pi_SC)
  mcw = (mcv - 1) / (sum(mcv) - length(mcv))   
  SSNV_model[["mut_class_w"]][["SM"]] = mcw[1]
  SSNV_model[["mut_class_w"]][["SC"]] = mcw[2]
##



   SSNV_model[["SSNV_skew"]] = SSNV_skew
# linear model of log(rho) = m*f_skew/2 + b
# fit from normal samples with out.p = 0.001
   b = -33.67232
   m = 82.77464  

# fit with out.p = 0.005
#   b = -45.86149
#   m = 110.51106 

   SSNV_model[["rho"]] = exp(m * (SSNV_model[["SSNV_skew"]]/2) + b)

   return(SSNV_model)
}


update_SSNV_mixture_weights = function( SSNV_model, som_mut_Q_tab, post_Prs ) 
{
  pi_som_theta_q = SSNV_model[["kPiSomThetaQ"]]

  mut_mat =  som_mut_Q_tab * matrix(post_Prs[,"Pr_somatic_clonal"], 
                                             nrow=nrow(som_mut_Q_tab), 
                                             ncol=ncol(som_mut_Q_tab), 
                                             byrow=TRUE) 

  nix = post_Prs[, "Pr_somatic_clonal"] == 0
  mut_mat = mut_mat[!nix, , drop=FALSE]
  if (nrow(mut_mat) > 0) {
    som_q = colSums(mut_mat, na.rm=TRUE)
  } else { 
    som_q = rep(0, length(pi_som_theta_q)) 
  }
  
## update Q multiplicity Dir mode
  som_theta_q_map = (pi_som_theta_q + som_q - 1) / (sum(pi_som_theta_q + som_q) - length(pi_som_theta_q)) 
  SSNV_model[["som_theta_Q_mode"]] = som_theta_q_map

## update mut_class_w Dir mode
  mut_w = post_Prs[, c("Pr_somatic_clonal", "Pr_subclonal"), drop=FALSE]
  mut_w = mut_w / rowSums(mut_w)
  mut_pr = colSums(mut_w, na.rm=TRUE)
  
  mut_class_w = SSNV_model[["mut_class_w"]]
  mcv = c(mut_class_w$Pi_SM, mut_class_w$Pi_SC)
  mcw = (mcv + mut_pr - 1) / (sum(mcv + mut_pr) - length(mcv))   
  mut_class_w$SM = mcw[1]
  mut_class_w$SC = mcw[2]
  SSNV_model[["mut_class_w"]] = mut_class_w

  return(SSNV_model)
}


H_exp_subclonal_SSNV = function( log.ccf.post, SSNV_model )
{
   ccf_grid = SSNV_model[["ccf_grid"]] 
   lambda = SSNV_model[["lambda"]]

   LL = dexp( ccf_grid, rate=lambda, log=TRUE) 
   LL = LL - LogAdd(LL)

   dx = diff( c(0,ccf_grid) )

   log.ev.H1 = grid_compose_LL( LL, NA, log.ccf.post, dx )

#   log.ev.H1 = LogAdd(post.numerator + log(dx))

   log.prior = log(SSNV_model[["mut_class_w"]][["SC"]] * SSNV_model[["W_E"]])

   return(log.ev.H1 + log.prior)
}

H_unif_subclonal_SSNV = function( log.ccf.post, SSNV_model )
{
   ccf_grid = SSNV_model[["ccf_grid"]] 

   LL = dunif( ccf_grid, 0, 1, log=TRUE) 
   LL = LL - LogAdd(LL)

   dx = diff( c(0,ccf_grid) )
   log.ev.H2 = grid_compose_LL( LL, NA, log.ccf.post, dx )

#   post.numerator = post.ccf.LL.grid + dunif( ccf_grid, 0, 1, log=TRUE) 
#   dx = diff(ccf_grid)
#   log.ev.H2 = LogAdd(post.numerator + log(dx))

   log.prior = log(SSNV_model[["mut_class_w"]][["SC"]] * SSNV_model[["W_U"]])

   return(log.ev.H2 + log.prior)
}




SSNV_FhatCombPost <- function(alt, ref, f.comb, SSNV_model) 
{
  cov = alt+ref
#  fhat= alt/cov

  f.comb = f.comb * SSNV_model[["SSNV_skew"]]
  log.pr <- matrix(NA, nrow = length(alt), ncol = length(f.comb))

  for (i in 1:length(f.comb)) 
  {
#    eqn. 12 Carter 2012
#    log.pr[, i] <- dbeta(f[i], fhat * cov + 1, (1 - fhat) * cov + 1, log = TRUE)
#    log.pr[, i] <- dbeta(f[i], s * fhat * cov + 1, (1 - s * fhat) * cov + 1, log = TRUE)

    f_skew = SSNV_model[["SSNV_skew"]]
    rho = SSNV_model[["rho"]]
    A = f_skew * f.comb[i] * rho
    B = (1-f_skew*f.comb[i]) * rho
    log.pr[,i] <- d_beta_binom( alt, A, B, cov, log=TRUE )
  }

#  if(any(is.nan(log.pr))) { stop("NaN log.pr!") }

  log.pr[is.nan(log.pr)] = 0  

  return(log.pr)
}


## integrate over subclonal SCNA f_c
mut_qm_fc_grid_integral = function(alt, ref, int_mat, outer_prod, SSNV_model) 
{
  cov = alt+ref
  h_qm_ll_mat = matrix(NA, nrow=length(alt), ncol=nrow(outer_prod))

  f_skew = SSNV_model[["SSNV_skew"]]
  rho = SSNV_model[["rho"]]
  A = f_skew * outer_prod * rho
  B = (1 - f_skew * outer_prod) * rho
  B[B < 0] = 0.0    ## fix rounding errors leading to small negative numbers

  for (j in 1:length(alt))
  {
    ldbb = d_beta_binom( alt[j], A, B, cov[j], log=TRUE )
#    ldbb[is.nan(ldbb)] = -Inf   ## can happen if alt and A are both 0
    ldbb[ outer_prod == 0 ] = ifelse( alt[j]==0, log(1), log(0) )
    ldbb[ outer_prod == 1 ] = ifelse( ref[j]==0, log(1), log(0) )

#    mat_sum = int_mat + ldbb
#    mat_sum[ int_mat== -Inf & ldbb == -Inf] = -Inf
#    h_qm_ll_mat[j, ] = LogAdd(mat_sum)
    h_qm_ll_mat[j, ] = LogAdd(int_mat + ldbb)

#   dbeta(outer_prod, fhat[j] * cov[j] + 1, (1 - fhat[j]) * cov[j] + 1, log=TRUE))
  }

  if( any(is.nan(h_qm_ll_mat))) { stop() }
  
  return(h_qm_ll_mat)
}

## double integral over subclonal SCNA f3 and f2 to get Pr(f2+f3)
joint_SSNV_SCNA_grid_density = function(alt, ref, AF_mat, joint_f2_f3_log_prior_mat, SSNV_model) 
{
  cov = alt+ref
  N_grid = length(SSNV_model[["ccf_grid"]])

#  SSNV_ccf_dens = matrix(NA, nrow=length(alt), ncol=N_grid )
#  H3_SSNV_ccf_dens = matrix(NA, nrow=length(alt), ncol=N_grid )
  f2_f3_Z = rep(NA, length(alt))
  f2_f3_dens = array( NA, dim=c(length(alt), N_grid, N_grid ) )

  f_skew = SSNV_model[["SSNV_skew"]]
  rho = SSNV_model[["rho"]]
  A = f_skew * AF_mat * rho
  B = (1 - f_skew * AF_mat) * rho
  B[B < 0] = 0.0    ## fix rounding errors leading to small negative numbers

  for (j in 1:length(alt))
  {
    ldbb = d_beta_binom( alt[j], A, B, cov[j], log=TRUE )

    ldbb[ AF_mat == 0 ] = ifelse( alt[j]==0, log(1), log(0) )
    ldbb[ AF_mat == 1 ] = ifelse( ref[j]==0, log(1), log(0) )
    ldbb[ AF_mat > 1 ] = -Inf

#    if(any(is.nan(ldbb))) { stop() }

    joint_LL =  joint_f2_f3_log_prior_mat + ldbb
    joint_LL[is.nan(joint_LL)] = -Inf  ## can happen because -Inf + -Inf = NaN

    f2_f3_Z[j] = LogAdd(as.vector(joint_LL))
    f2_f3_dens[j,,] = exp(joint_LL - f2_f3_Z[j])
  }

  if( any(is.nan(f2_f3_dens))) { stop() }

  return( list( "SSNV_log_Z"=f2_f3_Z, "f2_f3_dens"=f2_f3_dens) )
}


## Is this used?
calc_ccf_posterior_dens_grid = function(alt, ref, SSNV_model, alpha, q ) 
{
   ll_grid = calc_ccf_posterior_LL_grid(alt, ref, SSNV_model, alpha, q ) 
   pr_grid = exp(ll_grid - LogAdd(ll_grid))
   return(pr_grid)
}   

calc_ccf_posterior_LL_grid = function(alt, ref, SSNV_model, alpha, q ) 
{
## TODO:  add X-chr support
  f_skew =  SSNV_model[["SSNV_skew"]]
  rho =  SSNV_model[["rho"]]
  grid = SSNV_model[["ccf_grid"]]

  f = (alpha * grid) / (2 * (1 - alpha) + alpha * q)
  cov = alt + ref
#  ll_grid = dbinom(alt, cov, f, log=TRUE )

  use_f = f * f_skew
  A = use_f * rho
  B = (1-use_f) * rho

  ll_grid <- d_beta_binom( alt, A, B, cov, log=TRUE )

  if( length(alt)>1 )
  {
     ll_grid[ alt==0, f == 0 ] = log(1)
     ll_grid[ alt>0, f == 0 ] = log(0)

     ll_grid[ ref==0, f == 1 ] = log(1)
     ll_grid[ ref>0, f == 1 ] = log(0)
   }
   else
   {
     ll_grid[ f==0 ] = ifelse( alt==0, log(1), log(0) )
     ll_grid[ f==1 ] = ifelse( ref==0, log(1), log(0) )
   }

#  ll_grid[is.nan(ll_grid)] = -Inf   ## can happen if alt and A are both 0
  
  return(ll_grid)   
}


## TODO- reimplement for new model
CrypticScnaPost <- function(f.hat.vals, sDelta, cov, SSNV_model) {
  min.sc.af <- sDelta

  loglik <- rep(NA, length(f.hat.vals))

  beta.int <- pbeta(min.sc.af, f.hat.vals * cov + 1, (1 - f.hat.vals) * cov + 1,
                    lower.tail=FALSE, log.p=TRUE)
  loglik <- beta.int + log(1 / min.sc.af)
  
  return(loglik)
}

clonal_SSNV_Mult_Prior <- function(W, N, som.q.theta) {
  res <- som.q.theta[c(1:N)]
  res <- res / sum(res)
  return(W[["SM"]] * res)
}


GermlineMutPrior <- function(W) {
  res <- c(49, 49, 1)
  res <- res / sum(res)
  return(W[["GL"]] * res)
}


