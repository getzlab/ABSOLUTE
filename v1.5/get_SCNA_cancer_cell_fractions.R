get_subclonal_copy_states = function( obs, b, delta, SCNA_model, neg.ix )
{ 
  Q = SCNA_model[["kQ"]]

  comb_A =  GetCopyRatioComb(Q, delta, b, obs$error.model)
  comb_X =  GetCopyRatioComb(Q, delta, b/2, obs$error.model)
  
  ## convert dist over CR to dist over CCF for subclonal SCNAs - not dependant on theta.Q
  seg_q_tab = SCNA_model_clonal_seg_Q_post( obs, b, delta, SCNA_model )
#  modal_cn = which.max(colSums(seg_q_tab)) - 1
#  if(modal_cn==0) { modal_cn = 1 } ## don't model state of 0, screws up lots of calculations

  modal_cn = which.max(colSums(seg_q_tab)[-1] )  # find modal non-0 CN state

## Autosomes
  cr_set = obs$d.tx[!obs[["male_X"]] ] 
  res = calc_firstorder_subclonal_states(cr_set, modal_cn, comb_A, neg.ix )
  qc_A = res$qc
  qs_A = res$qs

## X-chr male?
  if( any(obs[["male_X"]]))
  {
    cr_set = obs$d.tx[obs[["male_X"]] ] 
    res = calc_firstorder_subclonal_states(cr_set, modal_cn, comb_X, neg.ix )
    qc_X = res$qc
    qs_X = res$qs

    qc = rep(NA, length(qc_A) + length(qc_X) )
    qs = rep(NA, length(qc_A) + length(qc_X) )

    qc[ !obs[["male_X"]] ] = qc_A
    qc[  obs[["male_X"]] ] = qc_X

    qs[ !obs[["male_X"]] ] = qs_A
    qs[  obs[["male_X"]] ] = qs_X
  }
  else {
    qc = qc_A
    qs = qs_A
  }

  return( cbind(qc, qs))
}



calc_SCNA_CCF_dens = function( obs, CN_states, b, delta, SCNA_model )
{
  Q = SCNA_model[["kQ"]]
  ccf_grid = SCNA_model[["ccf_grid"]]
  n_seg = length(obs[["d.tx"]])

  log_ccf_dens = matrix(NA, nrow=n_seg, ncol=length(ccf_grid))
  colnames(log_ccf_dens) = ccf_grid
  
  cr_dens = matrix(NA, nrow=n_seg, ncol=length(ccf_grid))
  cr_grid = matrix(NA, nrow=n_seg, ncol=length(ccf_grid))
  
  comb_A = GetCopyRatioComb(Q, delta, b, obs$error.model)
  comb_X = GetCopyRatioComb(Q, delta, b/2, obs$error.model)

  qc = CN_states[,"qc"]
  qs = CN_states[,"qs"]

  for (i in 1:nrow(log_ccf_dens)) 
  {
 ## outside comb domain: cr < comb(0) | cr > comb(Q) 
    if (is.na(qs[i])) { 
#      log_ccf_dens[i,] = log(SCNA_model[["Pr_sub_zero_CN"]])
      log_ccf_dens[i,] = NA
      next 
    }
    
    dd = (comb_A[qs[i] + 1] - comb_A[qc[i] + 1])
    cr_grid[i, ] = (dd * ccf_grid) + comb_A[qc[i] +1 ]
    
    cr_dens[i, ] = GetScnaStderrGridDensity(obs, cr_grid[i, ], SCNA_model, i )    
    log_ccf_dens[i, ] = cr_dens[i, ] - LogAdd(cr_dens[i, ]) 

    if(any( is.na(log_ccf_dens[i,]))) { stop("NA log_ccf_dens") }
  }

  return( log_ccf_dens )
}


# mode and 95-CI
calc_ccf_summary = function(ccf_dens)
{
  ccf_summary = matrix( NA, ncol=3, nrow=nrow(ccf_dens))

 ## filter out segs with any NAs in CCF dist
  na.ix = apply( is.na(ccf_dens), 1, any )
  ccf_dens = ccf_dens[!na.ix,, drop=FALSE]

  n_seg = nrow(ccf_dens)
  ccf_hat = rep(NA, n_seg)
  ccf_ci95 = matrix(NA, nrow=n_seg, ncol=2)

  ccf_grid = as.numeric(colnames(ccf_dens))

  for( i in 1:n_seg )
  {
    ccf_hat[i] = ccf_grid[which.max(ccf_dens[i, ])]
    ecdf = cumsum(ccf_dens[i, ])
    
    if(ecdf[1]==1) { ccf_ci95[i,]=c(0,0.1) }
    else
    {
       ccf_ci95[i, ] = approx(ecdf, y=ccf_grid, xout=c(0.025, 0.975))$y
    }
  }

  nix1 = is.na(ccf_ci95[, 1])
  ccf_ci95[nix1, 1] = min(ccf_grid)
  nix2 = is.na(ccf_ci95[, 2])
  ccf_ci95[nix2, 2] = max(ccf_grid)
  
  ## Round up in last bin.   TODO round down in 1st bin 
  ix = ccf_ci95[, 2] > ccf_grid[length(ccf_grid) - 1]
  ccf_ci95[ix, 2] = 1.0
  colnames(ccf_ci95) = c("CI95_low", "CI95_high")

  res=cbind( CCF_hat=ccf_hat, ccf_ci95 )
  ccf_summary[ !na.ix, ] = res

  return(ccf_summary)
}


calc_firstorder_subclonal_states = function(cr_set, modal_cn, comb, neg.ix ) 
{
  qs = rep(NA, length(cr_set))
  qc = rep(NA, length(cr_set))
  
  del_ix  = cr_set < comb[modal_cn+1]
  gain_ix = cr_set >= comb[modal_cn+1]
  Q = length(comb)
  
  ## qs = CN in the subclone with the SCNA
  del_vals = (c(0:(modal_cn - 1)))
  for (i in seq_along(del_vals)) {
    qs[del_ix & (cr_set > comb[del_vals[i] + 1]) ] = del_vals[i]
  }
  
 ## dels significantly below 0
  qs[del_ix & neg.ix] = NA    

 ## dels below 0 but clonal
  below0.ix = cr_set <= comb[1]
  qs[del_ix & !neg.ix & below0.ix ] = 0

  qc[ del_ix & !is.na(qs) ] = qs[del_ix & !is.na(qs)] + 1
  qc[ del_ix & is.na(qs)  ] = 0
  
  gain_vals = rev(c((modal_cn + 1) : Q ))
  for (i in seq_along(gain_vals)) {
    qs[gain_ix & (cr_set < comb[gain_vals[i] + 1])] = gain_vals[i]
  }

# amps above max comb have qs = NA
  qc[gain_ix & !is.na(qs)] = qs[gain_ix & !is.na(qs)] - 1
  qc[gain_ix &  is.na(qs)] = Q

  if( any(is.na(qc)) ) { stop() }
  
  if( any(qs > Q, na.rm=TRUE) ) { stop() }

  res = list(qs=qs, qc=qc)
  return(res)
}


## Not used currently
GetScnaPowTab <- function(obs, b, delta, CN_states, SCNA_model, alpha=0.05 ) 
{
   N <- length(obs[["d.tx"]])
   ccf_grid = SCNA_model[["ccf_grid"]]

   qc = CN_states[,"qc"]
   qs = CN_states[,"qs"]

    Q = SCNA_model[["kQ"]]
#    male_X = obs[["male_X"]]
    comb_A <- GetCopyRatioComb(Q, delta, b, obs[["error.model"]])

   seg_sigma = get_seg_sigma( SCNA_model, obs )

   pow = matrix(NA, nrow=N, ncol=length(ccf_grid))
   cr_grid = matrix(NA, nrow=N, ncol=length(ccf_grid))

   for( i in 1:N)
   {
      if(is.na(qc[i])) { pow[i,]=1; next}  ## sub-zero segs

      dd = (comb_A[qs[i] + 1] - comb_A[qc[i] + 1])
      cr_grid[i, ] = (dd * ccf_grid) + comb_A[qc[i] +1 ]

      zero = comb_A[qc[i]+1] ## dCR = 0

      y_05 = qnorm( 1-alpha, zero, seg_sigma[i])
      pow[i,] = 1 - pnorm( y_05, cr_grid[i,], seg_sigma[i] )
   } 

   return(pow)
}




GetScnaStderrGridDensity <- function(obs, grid, SCNA_model, i=NA) {
   N <- length(obs[["d.tx"]])
   grid.dens <- matrix(NA, nrow=N, ncol=length(grid))
   dx <- c(0, diff(grid))

   seg_sigma = get_seg_sigma( SCNA_model, obs )

   if( is.na(i))
   {
      for (i in seq_len(N)) {
        grid.dens[i, ] <- dnorm( grid, obs[["d.tx"]][i], seg_sigma[i], log=TRUE)
      }
   } else {
      grid.dens = dnorm( grid, obs$d.tx[i], seg_sigma[i], log=TRUE)
   }

   return(grid.dens)
}




reconcile_clonal_homdels_with_obs_SSNVs = function( mut.cn.dat, subclonal_ix, seg.q.tab, seg.post.subclonal, SCNA_model )
{
   clonal.mut.tab = get_muts_nearest_clonal_scna(mut.cn.dat, seg.q.tab, SCNA_model[["kQ"]])
## don't allow clonal homozygous deletions over regions with SSNVs having > 0 alt reads.
   homdel.ix = clonal.mut.tab[,"HS_q_hat_1"] == 0  &  clonal.mut.tab[,"HS_q_hat_2"] == 0 & mut.cn.dat[,"alt"] > 0

   if( any(homdel.ix) )
   {
      A1.ix = mut.cn.dat[homdel.ix,"A1.ix"]
      A2.ix = mut.cn.dat[homdel.ix,"A2.ix"]

     ## group to unique pairs of homologous segs:
      keys = unique( paste( A1.ix, A2.ix, sep="__") )

      for( i in 1:length(keys) )
      {
         res = strsplit(keys[i], "__" )[[1]]
         A1.ix = as.integer(res[1])
         A2.ix = as.integer(res[2])

         if( seg.post.subclonal[A1.ix] < seg.post.subclonal[A2.ix] )
         {
            subclonal_ix[A2.ix] = TRUE 
         }
         else
         {
            subclonal_ix[A1.ix] = TRUE 
         }
      }
   }
   return(subclonal_ix)
}

## Create SCNA clonality summary for SSNV models
allelic_get_subclonal_SCNA_info = function(  obs, b, delta, SCNA_model, mut.cn.dat )
{
   seg.post_mix.w = SCNA_model[["seg.post_mix.w"]]
   seg.q.tab = SCNA_model[["seg.q.tab"]]

   seg.post.subclonal = seg.post_mix.w[, "unif"]
   subclonal_ix = seg.post.subclonal > 0.1   ## used to select model for SSNVs


## For now, revert to unclustered CCF.   Consider using CCF from collapsed DP post ...
   log_ccf_dens = log( SCNA_model[["seg_CCF_dens"]] )
   ccf_sum = calc_ccf_summary( exp(log_ccf_dens) )


   if( !is.na(mut.cn.dat) )
   {
      subclonal_ix = reconcile_clonal_homdels_with_obs_SSNVs( mut.cn.dat, subclonal_ix, seg.q.tab, seg.post.subclonal, SCNA_model )
   }

   CN_states = SCNA_model[["CN_states"]]
   subclonal_scna_tab = cbind(ccf_sum, subclonal_ix, Pr_subclonal = seg.post.subclonal,  
				qs=CN_states[,"qs"], qc=CN_states[,"qc"] )

   return( list("subclonal_scna_tab"=subclonal_scna_tab, "log_ccf_dens"=log_ccf_dens) )
}


## Create SCNA clonality summary for SSNV models
total_get_subclonal_SCNA_info = function(  obs, b, delta, SCNA_model, mut.cn.dat )
{
   seg.post_mix.w = SCNA_model[["seg.post_mix.w"]]
   seg.q.tab = SCNA_model[["seg.q.tab"]]

   seg.post.subclonal = rowSums(seg.post_mix.w[, c("unif", "exp")]) 
   subclonal_ix = seg.post.subclonal > 0.1    ## used to select model for SSNVs

   CN_states = get_subclonal_copy_states( obs, b, delta, SCNA_model )
   log_ccf_dens = calc_SCNA_CCF_dens( obs, CN_states, b, delta, SCNA_model  )
   ccf_sum = calc_ccf_summary( exp(log_ccf_dens) )

   if( !is.na(mut.cn.dat))
   {
      clonal.mut.tab = get_muts_nearest_clonal_scna(mut.cn.dat, seg.q.hat, SCNA_model[["kQ"]])
## don't allow clonal homozygous deletions over regions with SSNVs having > 0 alt reads.
      homdel.ix = clonal.mut.tab[,"q_hat"] == 0  & mut.cn.dat[,"alt"] > 0

      if( any(homdel.ix) )
      {
         mut_seg_ix = mut.cn.dat[homdel.ix,"mut_seg_ix"]
         subclonal_ix[ mut_seg_ix ] = TRUE
      }
   }

   subclonal_scna_tab = cbind(ccf_sum, subclonal_ix, Pr_subclonal = seg.post.subclonal,  
				qs=CN_states[,"qs"], qc=CN_states[,"qc"] )

   return( list("subclonal_scna_tab"=subclonal_scna_tab, "log_ccf_dens"=log_ccf_dens) )
}



