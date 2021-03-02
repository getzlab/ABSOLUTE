SCNA_model_loglik <- function(obs, b, delta, SCNA_model, Provisional=FALSE )
{
  SCNA_seg_mix_LL = calc_SCNA_seg_mix_LL( obs, b, delta, SCNA_model, Provisional )
  SCNA_seg_LL = LogAdd(SCNA_seg_mix_LL)
  SCNA_LL = sum(SCNA_seg_LL) + SCNA_model_eval_prior(SCNA_model)

  if (!is.finite(SCNA_LL)) { stop(" Non-finite log-liklihood!") }

  return(SCNA_LL)
}

calc_SCNA_seg_mix_LL = function( obs, b, delta, SCNA_model, Provisional=FALSE, allelic_to_total_Wq0=FALSE )
{
 ## Clonal comb mixture
  log.clonal.seg.mat = SCNA_model_calc_log.clonal.seg.mat( obs, b, delta, SCNA_model, allelic_to_total_Wq0 )
  clonal_LL <- LogAdd(log.clonal.seg.mat)

 ## Amplified segs
  amp_LL <- SCNA_model_calc_Amp_LL( obs, b, delta, SCNA_model )
 

  log_mix_w = matrix(log(SCNA_model[["mix.w"]]),  nrow=nrow(log.clonal.seg.mat), ncol=length(SCNA_model[["mix.w"]]), byrow=TRUE )
  colnames(log_mix_w) = names(SCNA_model[["mix.w"]])

  if( !Provisional )
  {
     SC_LL_mat <- SCNA_model_calc_log.subclonal.seg.LL( obs, b, delta, SCNA_model )
  }
  else
  {
     SC_LL_mat = matrix( log(c(1e-2, 1e-2)), ncol=2, nrow=length(clonal_LL))
  }

  SCNA_seg_mix_LL <- cbind(clonal_LL, amp_LL, SC_LL_mat) + log_mix_w
  if( any( is.nan(SCNA_seg_mix_LL) ) ) {stop()}

  return(SCNA_seg_mix_LL)
}


SCNA_model_seg_mix.w_post = function( obs, b, delta, SCNA_model, allelic_to_total_Wq0=FALSE )
{
  SCNA_seg_mix_LL = calc_SCNA_seg_mix_LL( obs, b, delta, SCNA_model, allelic_to_total_Wq0 )
  seg_mix.w_post = exp(SCNA_seg_mix_LL - LogAdd(SCNA_seg_mix_LL))
  colnames(seg_mix.w_post) = names( SCNA_model[["pi.mix.w"]])

  return(seg_mix.w_post)
}

SCNA_model_clonal_seg_Q_post = function(  obs, b, delta, SCNA_model, allelic_to_total_Wq0=FALSE  )
{
  log.clonal.seg.mat = SCNA_model_calc_log.clonal.seg.mat( obs, b, delta, SCNA_model, allelic_to_total_Wq0 )
  seg_Q_post <- exp( log.clonal.seg.mat - LogAdd(log.clonal.seg.mat) )
  
  return(seg_Q_post)
}



# integrated composition of unif and exp over posterior seg CCF dist
SCNA_model_calc_log.subclonal.seg.LL = function( obs, b, delta, SCNA_model )
{
  dx = diff( c(0,SCNA_model[["ccf_grid"]]) )
  CN_states = get_subclonal_copy_states( obs, b, delta, SCNA_model )
  log.scna.ccf.post = calc_SCNA_CCF_dens( obs, CN_states, b, delta, SCNA_model )

  ## power for discrimination of obs CR from clonal at sig level alpha
#  log_SC_seg_pow_tab <- log(GetScnaPowTab(obs, b, delta, CN_states, SCNA_model, alpha=1e-2 ))
  log_SC_seg_pow_tab <- NA

# Uniform component
  dom = SCNA_model[["unif_CCF_dom"]]
  unif.LL <- dunif( SCNA_model[["ccf_grid"]], dom[1], dom[2], log=TRUE )
  unif.LL = unif.LL - LogAdd(unif.LL)
  subclonal_Unif_LL = grid_compose_LL( unif.LL, log_SC_seg_pow_tab, log.scna.ccf.post, dx )

# Exponential component
  exp.LL <- dexp( SCNA_model[["ccf_grid"]], rate=exp(SCNA_model[["Theta"]]["lambda"]), log=TRUE ) 
  exp.LL = exp.LL - LogAdd(exp.LL)
  subclonal_Exp_LL = grid_compose_LL( exp.LL, log_SC_seg_pow_tab, log.scna.ccf.post, dx )

  LL_seg_mat = cbind(subclonal_Exp_LL, subclonal_Unif_LL)

  return(LL_seg_mat)
}


SCNA_model_calc_log.clonal.seg.mat = function( obs, b, delta, SCNA_model, allelic_to_total_Wq0 )
{
  clonal_seg_LL = get_clonal_seg_LL_mat( obs, b, delta, SCNA_model )

## TODO: X-chr support!
  comb_A <- GetCopyRatioComb(SCNA_model[["kQ"]], delta, b, obs[["error.model"]])
  Wq0 = get_comb_Wq0( obs$e.cr, comb_A, SCNA_model[["Theta"]]["theta.Q"], SCNA_model[["theta.0"]])

  if( allelic_to_total_Wq0 == TRUE )  
  {
    ## create a prior for total CN states from the allelic prior
    tot.Wq0 = convolve( Wq0, rev(Wq0), type="open" )
    ## Truncate and renormalize
    tot.Wq0 = tot.Wq0[ c(1:SCNA_model[["kQ"]]) ]
    tot.Wq0[tot.Wq0 < .Machine$double.eps] = .Machine$double.eps  ## catch rounding errors from convolution
    tot.Wq0 = tot.Wq0 / sum(tot.Wq0)
    ## prior for CN=0 is not equivalent to convolution of allelic CN  (non-independant)  
    tot.Wq0[ c(2:SCNA_model[["kQ"]]) ] = tot.Wq0[ c(2:SCNA_model[["kQ"]]) ] / sum(tot.Wq0[ c(2:SCNA_model[["kQ"]]) ]) * (1- 1e-3) 
    tot.Wq0[1] = 1e-3  ## Hacky way of making comb prior consistent with total-CN-only version of model

    Wq0 = tot.Wq0
  }

  log.seg.pQ = matrix(log(Wq0), nrow=length(obs$d.tx), ncol=SCNA_model[["kQ"]], byrow=TRUE )
  log.clonal.seg.mat <- clonal_seg_LL + log.seg.pQ

  if( any(!is.finite(log.clonal.seg.mat))){stop()}

  return(log.clonal.seg.mat)
}


SCNA_model_calc_Amp_LL = function( obs, b, delta, SCNA_model )
{
   res = segs_in_comb_domain( obs, b, delta, SCNA_model )
   amp_ix = res[["amp_ix"]]

   LL = rep( -Inf, length(obs[["d.tx"]]))
   LL[amp_ix] = 0

   return(LL)
}
 


get_clonal_seg_LL_mat = function( obs, b, delta, SCNA_model )
{
  d = obs$d.tx
  use_sigma = get_seg_sigma( SCNA_model, obs )

  Q = SCNA_model[["kQ"]]
  male_X = obs[["male_X"]]

  comb_A <- GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
  clonal_seg_LL_A = sapply(comb_A, dnorm, d[!male_X], use_sigma[!male_X], log=TRUE)

  if( any(male_X))
  {
    comb_X <- GetCopyRatioComb(Q, delta, b/2, obs[["error.model"]])
    clonal_seg_LL_X = sapply(comb_X, dnorm, d[male_X], use_sigma[male_X], log=TRUE)

    clonal_seg_LL = matrix( NA, nrow=length(d), ncol=Q )
    clonal_seg_LL[ male_X, ] = clonal_seg_LL_X
    clonal_seg_LL[ !male_X, ] = clonal_seg_LL_A
  }
  else{ clonal_seg_LL = clonal_seg_LL_A }

  if(any(is.na(clonal_seg_LL))) { stop() }

  return(clonal_seg_LL)
}


get_comb_interval_gaussian_CDF_seg_mat = function( obs, b, delta, SCNA_model )
{
  get_int_mat = function( d, comb, sigma )
  {
     Q = length(comb)
     seg_interval_dens_mat = matrix( NA, nrow=length(d), ncol=Q+1 )

     for( i in 1:length(d))
     {
        mid = diff( pnorm(comb, d[i], sigma[i] ) )
        neg = pnorm( comb[1], d[i], sigma[i] )
        amp = 1 - pnorm( comb[Q], d[i], sigma[i] )

        seg_interval_dens_mat[i,] = c( neg, mid, amp )
     }

     return( seg_interval_dens_mat )
  }

  d = obs[["d.tx"]]
  use_sigma = get_seg_sigma( SCNA_model, obs )
  Q = SCNA_model[["kQ"]]
  male_X = obs[["male_X"]]
  comb_A <- GetCopyRatioComb(Q, delta, b, obs[["error.model"]])

  seg_interval_dens_mat = matrix( NA, nrow=length(d), ncol=Q+1 )

  seg_interval_dens_mat[ !male_X, ] = get_int_mat( d[!male_X], comb_A, use_sigma[!male_X] )

  if( any(male_X)) 
  {
     comb_X <- GetCopyRatioComb(Q, delta, b/2, obs[["error.model"]])
     seg_interval_dens_mat[ male_X, ] = get_int_mat( d[male_X], comb_X, use_sigma[male_X] )  
  }

  return(seg_interval_dens_mat)
}



## Gaussian prior on CN states centered at e.cr, with variance theta.Q and comb[0] = theta.0
get_comb_Wq0 = function( e.cr, comb, theta.Q, theta.0 )
{
   if(max(comb)<1){ use_out = max(comb) }
   else{ use_out = 1 }

   mu = approx( y=c(0:(length(comb)-1)), x=comb, xout=use_out)$y
   Q = length(comb)
   sigma = exp(theta.Q)
   pq = dnorm( c(1:(Q-1)), mu, sigma, log=FALSE )

   if( !is.na(theta.0)) 
   {
      p0 = theta.0
   }
   else
   {
      p0 = dnorm( 0, mu, sigma, log=FALSE )
   }

   epsilon = 1e-10
   Wq0 = c( p0, pq ) + epsilon

   if( any(!is.finite(Wq0))) { stop() }

   Wq0 = Wq0 / sum(Wq0)

   return( Wq0 )
}


