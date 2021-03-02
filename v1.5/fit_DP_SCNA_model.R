SCNA_model_calc_CCF_DP_loglik = function( obs, b, delta, SCNA_model, verbose=verbose )
{
   comb_A <- GetCopyRatioComb(SCNA_model[["kQ"]], delta, b, obs[["error.model"]])

# Annotate segs for comb domain
   out_Prs = get_neg_and_amp_Prs( obs, b, delta, SCNA_model )
   neg_Pr = out_Prs[,"neg_Pr"]
   SCNA_model[["neg.Pr"]] = neg_Pr
   neg.ix = SCNA_model[["neg.ix"]] = neg_Pr > 0.99   #(1 - SCNA_model[["Pr_sub_zero_CN"]])

## Estimate derived/ancestor CN states and compute SCNA CCF distributions
   CN_states = get_subclonal_copy_states( obs, b, delta, SCNA_model, neg.ix )
   SCNA_model[["CN_states"]] = CN_states
   SCNA_model[["seg_CCF_dens"]] = exp( calc_SCNA_CCF_dens( obs, CN_states, b, delta, SCNA_model  ) )

   amp.ix = segs_in_comb_domain( obs, b, delta, SCNA_model )[["amp_ix"]]

   Pr_clonal = rowSums( SCNA_model[["seg_CCF_dens"]][ , SCNA_model[["clonal_CCF_bins"]] ] )
   clonal.ix = Pr_clonal > SCNA_model[["provisional_clonal_Pr_threshold"]] & !(amp.ix | neg.ix)  # provisional threshold
   SCNA_model[["bi.W.ix"]] = obs[["bi.allelic"]] # & obs[["W"]] < 1e-2
   high.sem.ix = obs[["d.stderr"]] > SCNA_model[["seg_sem_thresh"]] & !(amp.ix | neg.ix | clonal.ix )
# TODO - remove unpowered segs instead
#  log_SC_seg_pow_tab <- log(GetScnaPowTab(obs, b, delta, CN_states, SCNA_model, alpha=1e-2 ))

   SCNA_model[["seg.ix.tab"]] = cbind(amp.ix, neg.ix, clonal.ix, high.sem.ix, "bi.W.ix"=SCNA_model[["bi.W.ix"]] )

   SCNA_model[["collapsed_CCF_dens"]] = cbind( Pr_clonal, SCNA_model[["seg_CCF_dens"]][, - SCNA_model[["clonal_CCF_bins"]] ] )
   nc = ncol(SCNA_model[["collapsed_CCF_dens"]])
   collapsed_dx = rep( 1/nc, nc )

# Kernel density estimator
   SCNA_model[["collapsed_CCF_KDE"]] = colMeans(SCNA_model[["collapsed_CCF_dens"]][!(amp.ix|neg.ix),,drop=FALSE])
# Run DP mcmc
   SCNA_model = run_seg_CCF_DP( SCNA_model )

## Compose seg CR gaussians and DP predictive density estimate over CCF on collapsed CCF density 
   DP.out.ix = SCNA_model[["DP.out.ix"]]

   seg_CCF_LL = rep(NA, length(obs[["d.tx"]]))
   seg_CCF_LL[!DP.out.ix] = grid_compose_LL( log(SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]]), NA, log(SCNA_model[["collapsed_CCF_dens"]][!DP.out.ix,,drop=FALSE]), collapsed_dx )
   seg_CCF_LL[high.sem.ix] = grid_compose_LL( log(SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]]), NA, log(SCNA_model[["collapsed_CCF_dens"]][high.sem.ix,, drop=FALSE]), collapsed_dx )


## CN component of seg LL
   CN.seg.ix = !(amp.ix | neg.ix)
   res = optimize_theta.Q( obs, comb_A, SCNA_model, CN_states, CN.seg.ix )
   SCNA_model[["Theta"]][["theta.Q"]] = res[["theta.Q"]]
   seg_CN_LL = res[["seg_CN_LL"]]

## For segs in comb domain, likelihood is product of Gaussian on derived integer CN (mu=e.cr, sigma=theta.Q), and CCF on empirical (DP) density estimate
   seg_LL = rep(NA, length(obs[["d.tx"]]))
   seg_LL[!DP.out.ix] = seg_CN_LL[!DP.out.ix] + seg_CCF_LL[!DP.out.ix]
  
  ## Segs in comb dom, but not used in DP.    
   seg_LL[ clonal.ix ] = seg_CN_LL[clonal.ix] + 
                         log(SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]][1])    #clonal bin of predictive density estimate

   seg_LL[ high.sem.ix ] = seg_CN_LL[high.sem.ix] + seg_CCF_LL[high.sem.ix]
 
## For segs outside comb domain, add terms for neg/amp segs.
## Hard-settings for segs out of domain - TODO - make LL dependent on SCNA seg size
## FOR now - handle small bi-allelic events specially - they are probably artifacts - 
## TODO: compromise 
   bi.W.ix = SCNA_model[["bi.W.ix"]] 

#   M_b = SCNA_model[["amp_seg_Pr_M_b"]]
#   amp_Pr = M_b[1] / obs[["W"]][amp.ix & !bi.W.ix] + M_b[2]
#   eps = 1e-10
#   amp_Pr = 1 - obs[["W"]][amp.ix & !bi.W.ix]^(1/2) + eps
   if( any( amp.ix & !bi.W.ix )) 
   {
      amp_Pr = dexp(obs[["W"]][amp.ix & !bi.W.ix], 500) / dexp(0,500) 
      seg_LL[ amp.ix & !bi.W.ix ] = log( amp_Pr )
   }

   neg.rate = 2500
   if( any( neg.ix & !bi.W.ix )) 
   {
      neg_Pr = dexp(obs[["W"]][neg.ix & !bi.W.ix], neg.rate) / dexp(0,neg.rate)
      seg_LL[ neg.ix & !bi.W.ix ] = log(neg_Pr)
   }

   seg_LL[ (amp.ix | neg.ix) & bi.W.ix ] = log(SCNA_model[["Pr_bi.allelic.outlier"]])
   if( any( is.na(seg_LL) ) ) { stop() }

   SCNA_model[["seg_LL"]] = seg_LL
   SCNA_model[["seg_CN_LL"]] = seg_CN_LL
   SCNA_model[["LL"]] =  sum(seg_LL) 

   print(paste("loglik = ", round(SCNA_model[["LL"]] ,4), sep=""))
   print_model(SCNA_model)

   SCNA_model = get_seg_mix_posts( obs, SCNA_model )

   return(SCNA_model)
}


skip_seg_CCF_DP = function( SCNA_model )
{
## Handle rare edge case where no segs are left for DP clustering (e.g. all clonal)
## Fufill run_seg_CCF_DP contract with SCNA_model without actually running the DP

   nc = ncol( SCNA_model[["collapsed_CCF_dens"]] )
   clonal_dens = rep(0, nc)
   clonal_dens[1] = 1

   SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]] = clonal_dens
   SCNA_model[["seg_CCF_DP"]][["predictive_density_summary"]] = matrix( clonal_dens, nrow=4, ncol=length(clonal_dens), byrow=TRUE )

   SCNA_model[["seg_CCF_DP"]][["use_k"]] = 1
   SCNA_model[["seg_CCF_DP"]][["tree_clust"]] = list( "assign"=NA, "CCF_dens"=matrix(clonal_dens, nrow=1))


## For segs excluded due to high sem, fill in a regularized CCF density estimate
   collapsed_DP_CCF_dens = matrix( NA, nrow=nrow(SCNA_model[["seg_CCF_dens"]]), ncol=nc)
   high.sem.ix =  SCNA_model[["seg.ix.tab"]][, "high.sem.ix"]
   CCF_LL = t( log(SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]]) + t(log(SCNA_model[["collapsed_CCF_dens"]][high.sem.ix,]) ) ) 
   CCF_dens = exp(CCF_LL-LogAdd(CCF_LL))
   collapsed_DP_CCF_dens[high.sem.ix,] = CCF_dens
   SCNA_model[["seg_CCF_DP"]][["collapsed_DP_CCF_dens"]] = collapsed_DP_CCF_dens


   return(SCNA_model)
}


run_seg_CCF_DP = function( SCNA_model )
{
   collapsed_CCF_dens = SCNA_model[["collapsed_CCF_dens"]]

# all exclusive reasons to remove segs from DP clustering
   DP.out.ix = apply( SCNA_model[["seg.ix.tab"]][ , c("amp.ix", "neg.ix", "clonal.ix", "high.sem.ix") ], 1, any )
   SCNA_model[["DP.out.ix"]] = DP.out.ix
   N_DAT = sum(!DP.out.ix)

   if( N_DAT <= 1 )  ## No segs for clustering!!
   {
      ## if only a single segment is left, make sure it is high.sem, so it gets taken care of without DP.
      if(N_DAT == 1) 
      {
         s.ix = which( !DP.out.ix )
         SCNA_model[["seg.ix.tab"]][s.ix, "high.sem.ix"] = TRUE
      }

      print( "Skipping DP clustering...")
      return( skip_seg_CCF_DP( SCNA_model ))
   }

## set up pseudo-counts for clonal segs
   clonal.ix =  SCNA_model[["seg.ix.tab"]][, "clonal.ix"]
   high.sem.ix =  SCNA_model[["seg.ix.tab"]][, "high.sem.ix"]
   npc = sum(clonal.ix)
   pseudo_obs = matrix( 0, nrow=npc, ncol=ncol(collapsed_CCF_dens))
   pseudo_obs[,1] = 1

## run gibbs sampler to fit DP. 
   DP_input_densities = collapsed_CCF_dens[!DP.out.ix,, drop=FALSE]
   nc = ncol(DP_input_densities)
# Create a uniform prior for CCF loc
   clust_CCF_prior = rep(1/nc, nc)

   set.seed(0)  ## reproducable results!
   if( nrow(pseudo_obs) > 0 )
   {
      fixed_dens = pseudo_obs[1,,drop=FALSE]
      fixed_scales = npc
      DP_res = run_DP( DP_input_densities, SCNA_model[["N_DP_iter"]], SCNA_model[["Pi_k"]], clust_CCF_prior, use_fixed=TRUE, fixed_dens=fixed_dens, fixed_scales=fixed_scales ) 
      pred = summarize_predictive_density( DP_res[["DP_post"]], SCNA_model[["N_DP_burn"]], c(1:nc), use_fixed=TRUE, fixed_dens=fixed_dens, fixed_scales=fixed_scales )
   }
   else
   {
      DP_res = run_DP( DP_input_densities, SCNA_model[["N_DP_iter"]], SCNA_model[["Pi_k"]], clust_CCF_prior, use_fixed=FALSE )
      pred = summarize_predictive_density( DP_res[["DP_post"]], SCNA_model[["N_DP_burn"]], c(1:nc), use_fixed=FALSE )
   }

   ## estimate predictive density over CCF
   SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]] = pred[4,]  ## mean
   SCNA_model[["seg_CCF_DP"]][["predictive_density_summary"]] = pred

  ## summarize number of DP components
   est_k = summarize_DP_num_components( DP_res, SCNA_model[["N_DP_burn"]] )
   res=  rle( sort(est_k))
   rr = cbind( res[["lengths"]], res[["values"]] )
   r.ix = order( res[["values"]] )
   rr = rr[r.ix, , drop=FALSE ]
   if(nrow(rr)==1) { use_k = rr[1,2] }
   else {
      if( rr[1,1] * 10 < rr[2,1] ) { use_k = rr[2,2] } else { use_k = rr[1,2] }
   }
   SCNA_model[["seg_CCF_DP"]][["use_k"]] = use_k

  ## co-phenetic clustering of DP gibbs sampler iterations
   if( nrow(pseudo_obs) > 0 )
   {
      tree_clust = tree_cluster_DP( DP_input_densities, DP_res, SCNA_model[["N_DP_burn"]], use_k, clust_CCF_prior,
                                    use_fixed=TRUE, fixed_dens=fixed_dens, fixed_scales=fixed_scales )
   }
   else
   {
      tree_clust = tree_cluster_DP( DP_input_densities, DP_res, SCNA_model[["N_DP_burn"]], use_k, clust_CCF_prior )
   }
   SCNA_model[["seg_CCF_DP"]][["tree_clust"]] = tree_clust

   if( any(is.na(tree_clust)) ) { stop("NA in tree_clust")}  

#   tree_clust = resort_tree_clusters( SCNA_model, tree_clust )
#   SCNA_model[["seg_CCF_DP"]][["tree_clust"]] = tree_clust
#   SCNA_model[["seg_CCF_DP"]][["seg_clust_tab"]] = get_seg_clust_tab( SCNA_model )

## Summarize segs CCF dist, post DP
   collapsed_DP_CCF_dens = matrix( NA, nrow=nrow(SCNA_model[["seg_CCF_dens"]]), ncol=nc)

   res = summarize_mut_locations( DP_res[["DP_post"]], SCNA_model[["N_DP_burn"]] )
   collapsed_DP_CCF_dens[!DP.out.ix,] = res

## For segs excluded due to high sem, fill in a regularized CCF density estimate
   CCF_LL = t( log(SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]]) + t(log(collapsed_CCF_dens[high.sem.ix,, drop=FALSE]) ) ) 
   CCF_dens = exp(CCF_LL-LogAdd(CCF_LL))
   collapsed_DP_CCF_dens[high.sem.ix,] = CCF_dens

   if( any( !is.finite(collapsed_DP_CCF_dens[!DP.out.ix | high.sem.ix,])) ) { stop() }
   SCNA_model[["seg_CCF_DP"]][["collapsed_DP_CCF_dens"]] = collapsed_DP_CCF_dens

   return(SCNA_model) 
}

resort_tree_clusters = function ( SCNA_model, tree_clust )
{
## resort clusters by CCF
   mode_CCF_ix = apply(tree_clust[["CCF_dens"]], 1, which.max )

   clonal.ix =  SCNA_model[["seg.ix.tab"]][, "clonal.ix"]
#   high.sem.ix =  SCNA_model[["seg.ix.tab"]][, "high.sem.ix"]
   N = length(clonal.ix)

   assign = rep(NA, length(clonal.ix) )
   nix = apply( SCNA_model[["seg.ix.tab"]][ , c("amp.ix", "neg.ix", "high.sem.ix") ], 1, any )
   assign[ !clonal.ix & !nix ] = tree_clust[["assign"]]

##
   o.ix = order(mode_CCF_ix)
   assign[ clonal.ix ] = o.ix[1]   ## assignment to clonal cluster

   nix = is.na(assign)
   new.assign = rep(NA, length(clonal.ix) )
   new.assign[!nix] = o.ix[ assign[!nix] ]

#   tree_clust[["CCF_dens"]] = tree_clust[["CCF_dens"]][o.ix,]
#   tree_clust[["assign"]] = new.assign
   tree_clust[["assign"]] = assign

   tree_clust[["CCF_order"]] = o.ix

   return( tree_clust )
}


get_seg_clust_tab = function( SCNA_model )
{
   tree_clust = SCNA_model[["seg_CCF_DP"]][["tree_clust"]]
   K = nrow( tree_clust[["CCF_dens"]] )
   N = nrow( SCNA_model[["seg.ix.tab"]] )

## get a table with distr over components for all segs included
   seg_clust_tab = matrix( NA, nrow=N, ncol=K)

   log_data = log(SCNA_model[["collapsed_CCF_dens"]])
   log_clust_dens = log( tree_clust[["CCF_dens"]] )
   n.ix = apply( SCNA_model[["seg.ix.tab"]][ , c("amp.ix", "neg.ix" ) ], 1, any )

   for( n in 1:N )
   {
      if( n.ix[n] ) { next }

      for( k in 1:K )
      {
         seg_clust_tab[n,k] = log_cluster_prob( log_clust_dens[k,, drop=FALSE], log_data, n )
      }
   }

   clust_mix_W = rep(NA, K)
   for( k in 1:K )
   {
      clust_mix_W[k] = sum( tree_clust[["assign"]] == k, na.rm=TRUE )
   }
   clust_mix_W = clust_mix_W / sum(clust_mix_W)

   seg_clust_tab = seg_clust_tab + matrix( log(clust_mix_W), ncol=K, nrow=N, byrow=TRUE )
   seg_clust_tab = exp( seg_clust_tab - LogAdd(seg_clust_tab))

# For debugging - check consistency with original assignments
   seg_modal_clust = rep(NA, nrow(seg_clust_tab))
   seg_modal_clust[!n.ix] = apply(seg_clust_tab[!n.ix, , drop=FALSE], 1, which.max)

#   same = sum(tree_clust[["assign"]][!n.ix] == seg_modal_clust, na.rm=T)
#   bix = which( seg_modal_clust != tree_clust[["assign"]] )

   return( seg_clust_tab ) 
}




optimize_theta.Q = function( obs, comb, SCNA_model, CN_states, use.seg.ix )
{
   calc_seg_CN_LL = function( obs, comb, theta.Q, use.seg.ix, SCNA_model, CN_states )
   {
     seg_CN_LL = rep(NA, length(obs[["d.tx"]]))

     Wq0 = get_comb_Wq0( obs$e.cr, comb, theta.Q, SCNA_model[["theta.0"]])

## use.seg.ix is TRUE for segs in comb domain...
## Assume all segs are in derived CN state unless NA
     use_state = rep(NA, length(use.seg.ix) )
     uCN =  CN_states[use.seg.ix,2]
     if(any(is.na(uCN))) { stop() }
     use_state[use.seg.ix] = uCN

#     use_state[is.na(use_state)] = CN_states[use.seg.ix,1][is.na(use_state)]

## put amps at end of comb
#     use_state[ use_state >= length(Wq0) ] = length(Wq0)-1

     seg_CN_LL[use.seg.ix] = log( Wq0[use_state[use.seg.ix]+1] )    ## Wq0[1] is for CN=0

     if( any( is.na( seg_CN_LL[use.seg.ix] ))) { stop() }

     return(seg_CN_LL)
   }

   
   LL_wrap <- function(par)
   {
     theta.Q = par

     seg_CN_LL = calc_seg_CN_LL( obs, comb, theta.Q, use.seg.ix, SCNA_model, CN_states )
     LL = sum( seg_CN_LL[use.seg.ix] )

     if (!is.finite(LL)) 
     {
         print(paste(parname, ": Non-finite log-liklihood!", sep=""))
         stop()
     }
     else
     {
        if (verbose) { cat(symbol) }
     }

     return(-LL)
   }


## Calls calc_seg_CN_LL( obs, comb_A, theta.Q, use.seg.ix, SCNA_model, CN_states ) to optimize theta.Q
   theta.ix = which(names(SCNA_model[["Theta"]])=="theta.Q")
   limits = list("lower"=SCNA_model[["Theta_lower"]][theta.ix], "upper"=SCNA_model[["Theta_upper"]][theta.ix])
   symbol = SCNA_model[["Theta_sym"]][theta.ix]
   theta.Q = optimize(LL_wrap, lower=limits[["lower"]], upper=limits[["upper"]], tol=SCNA_model[["opttol"]], maximum=FALSE )[["minimum"]]

## compute copy-state LL using fit theta.Q (hat)
   seg_CN_LL = calc_seg_CN_LL( obs, comb, theta.Q, use.seg.ix, SCNA_model, CN_states )

   return( list("theta.Q"=theta.Q, "seg_CN_LL"=seg_CN_LL) )
}


get_seg_mix_posts = function( obs, SCNA_model )
{
   clonal.ix = SCNA_model[["seg.ix.tab"]][ , "clonal.ix"] 
   amp.ix = SCNA_model[["seg.ix.tab"]][ , "amp.ix"] 
   neg.ix = SCNA_model[["seg.ix.tab"]][ , "neg.ix"] 

## assignment of segs to mixtures
   seg.post_mix.w = matrix( 0, nrow=length(obs[["d.tx"]]), ncol=4 )

   colnames(seg.post_mix.w) = names(SCNA_model[["pi.mix.w"]])
   seg.post_mix.w[ amp.ix, "amp" ] = 1
   seg.post_mix.w[ neg.ix, "unif" ] = 1
   seg.post_mix.w[, "exp" ] = 0

## for collapsed seg CCF model:
   seg.post_mix.w[ clonal.ix, "clonal"] = 1  ## presumed clonal prior to clustering
   subclone_eval_ix = !(clonal.ix | amp.ix | neg.ix)
   seg.post_mix.w[subclone_eval_ix, "clonal"] = SCNA_model[["seg_CCF_DP"]][["collapsed_DP_CCF_dens"]][subclone_eval_ix,1]
   seg.post_mix.w[ subclone_eval_ix, "unif"] = 1 - seg.post_mix.w[ subclone_eval_ix, "clonal"]
   seg.post_mix.w[ subclone_eval_ix, "unif"][ seg.post_mix.w[ subclone_eval_ix, "unif"] < 0 ] = 0   ## rounding error
   SCNA_model[["seg.post_mix.w"]] = seg.post_mix.w

#  "mix.w": Hierarchical mixture weights - used as a prior for each seg to weigh LL and compute subclonal, amp Prs
#  names(SCNA_model[["pi.mix.w"]]) =c("clonal", "amp", "exp", "unif")
#  exp is not used - only unif for subclonal weight
# Don't bother fitting mix.w as hierarchical param now
   SCNA_model[["mix.w"]] = colMeans( SCNA_model[["seg.post_mix.w"]] )
   names(SCNA_model[["mix.w"]]) = names(SCNA_model[["pi.mix.w"]])

   if( any(!is.finite(SCNA_model[["seg.post_mix.w"]]))) { stop() }

   return(SCNA_model)
}


get_neg_and_amp_Prs = function( obs, b, delta, SCNA_model )
{
   seg_interval_dens_mat = get_comb_interval_gaussian_CDF_seg_mat( obs, b, delta, SCNA_model )

   neg_Pr = seg_interval_dens_mat[, 1]
   amp_Pr = seg_interval_dens_mat[, ncol(seg_interval_dens_mat)]

   return( cbind(neg_Pr, amp_Pr) )
}



segs_in_comb_domain = function( obs, b, delta, SCNA_model )
{
   comb_A = GetCopyRatioComb(SCNA_model[["kQ"]], delta, b, obs[["error.model"]])
   comb_X = GetCopyRatioComb(SCNA_model[["kQ"]], delta, b/2, obs[["error.model"]])
   male_X = obs[["male_X"]]

   amp_ix = rep(NA, length(obs[["d.tx"]]))
   neg_ix = rep(NA, length(obs[["d.tx"]]))

   amp_ix[!male_X] = obs[["d.tx"]] > max(comb_A)
   neg_ix[!male_X] = obs[["d.tx"]] < min(comb_A)

   amp_ix[male_X] = obs[["d.tx"]] > max(comb_X)
   neg_ix[male_X] = obs[["d.tx"]] < min(comb_X)

   return( list("amp_ix"=amp_ix, "neg_ix"=neg_ix))
}


apply_allelic_model_to_total_copy_ratios = function( tot.obs, b, delta, SCNA_model )
{
#   hom.seg.pair.tab = tot.obs[["hom.seg.pair.tab"]]  ## mapping to allelic segs
   seg.post_mix.w = SCNA_model[["seg.post_mix.w"]] ## for allelic data

 ## Hack to override 'automatic' selection of correct comb function which is set to "allelic" function at this point
   GetCopyRatioComb <<- total_get_copy_ratio_comb
  ##
   SCNA_model[["tot.seg.q.tab"]] = SCNA_model_clonal_seg_Q_post( tot.obs, b, delta, SCNA_model, allelic_to_total_Wq0=TRUE )

   n_tot = nrow( SCNA_model[["tot.seg.q.tab"]] )   #nrow( seg_pair_tab )
   tot.seg.post_mix.w = matrix( 0, nrow=n_tot, ncol=ncol(seg.post_mix.w) )
   colnames(tot.seg.post_mix.w) = colnames(seg.post_mix.w)

# Annotate segs for comb domain
   out_Prs = get_neg_and_amp_Prs( tot.obs, b, delta, SCNA_model )
   neg_Pr = out_Prs[,"neg_Pr"]
   SCNA_model[["neg.Pr"]] = neg_Pr
   neg.ix = SCNA_model[["tot.neg.ix"]] = neg_Pr > 0.999 
   amp.ix = segs_in_comb_domain( tot.obs, b, delta, SCNA_model )[["amp_ix"]]
   tot.seg.post_mix.w[amp.ix, "amp"] = 1

## Estimate derived/ancestor CN states and compute SCNA CCF distributions
   tot.CN_states = get_subclonal_copy_states( tot.obs, b, delta, SCNA_model, neg.ix )
   seg_CCF_dens = exp( calc_SCNA_CCF_dens( tot.obs, tot.CN_states, b, delta, SCNA_model  ) )

## Calc subclonal seg LL using DP CCF model
   Pr_clonal = rowSums( seg_CCF_dens[ , SCNA_model[["clonal_CCF_bins"]] ] )
   clonal.ix = Pr_clonal > 0.1 & !(amp.ix | neg.ix)  # provisional threshold
   tot.seg.post_mix.w[ clonal.ix, "clonal"] = 1  ## presumed clonal prior to clustering

# regularize total CN CCF estimates using predictive DP density learned from allelic data
   tot_collapsed_CCF_dens = cbind( Pr_clonal, seg_CCF_dens[, - SCNA_model[["clonal_CCF_bins"]] ] )
   CCF_LL = log(SCNA_model[["seg_CCF_DP"]][["CCF_DP_dens_est"]]) + log(tot_collapsed_CCF_dens ) 
   CCF_dens = exp(CCF_LL-LogAdd(CCF_LL))
   SCNA_model[["tot_collapsed_DP_CCF_dens"]] = CCF_dens


   subclone_eval_ix = !(clonal.ix | amp.ix | neg.ix)
   tot.seg.post_mix.w[ subclone_eval_ix, "clonal"] = SCNA_model[["tot_collapsed_DP_CCF_dens"]][subclone_eval_ix,1]
   tot.seg.post_mix.w[ subclone_eval_ix, "unif"] = 1 - tot.seg.post_mix.w[ subclone_eval_ix, "clonal"]
   tot.seg.post_mix.w[ subclone_eval_ix, "unif"][ tot.seg.post_mix.w[ subclone_eval_ix, "unif"] < 0 ] = 0   ## rounding error
   SCNA_model[["tot.seg.post_mix.w"]] = tot.seg.post_mix.w

   SCNA_model[["tot.seg.z.tab"]] = rowSums(tot.seg.post_mix.w[, c("unif", "exp")]) 
   SCNA_model[["tot.seg.amp.tab"]] = tot.seg.post_mix.w[,"amp"] 

   fac = (1-SCNA_model[["tot.seg.z.tab"]])
   fac[ fac < 0] = 0  ## round error
   SCNA_model[["tot.seg.qz.tab"]] = cbind( fac*SCNA_model[["tot.seg.q.tab"]], SCNA_model[["tot.seg.z.tab"]] )
 
  ## set back to allelic version
   GetCopyRatioComb <<- AllelicGetCopyRatioComb


## TODO: Also compute from allelic version, evaluate consistency
#   na.ix = apply( is.na(hom.seg.pair.tab), 1, any ) ## filtered allelic segs

#   tot.seg.post_mix.w[!na.ix, "clnal"] = seg.post_mix.w[ seg_pair_tab[,1], "clonal" ]   *   seg.post_mix.w[ seg_pair_tab[,2], "clonal" ]
#   tot.seg.post_mix.w[!na.ix, "amp"] = seg.post_mix.w[ seg_pair_tab[,1], "amp" ] + seg.post_mix.w[ seg_pair_tab[,2], "amp" ]
#   tot.seg.post_mix.w[!na.ix, "unif"] = 1 - rowSums( tot.seg.post_mix.w )
 

   return( SCNA_model )
}



calc_mode_seg_tabs = function( seg.obj, SCNA_model, b, delta )
{
   obs <- seg.obj[["obs.scna"]]

#   SCNA_model = get_seg_mix_posts( obs, SCNA_model )

   SCNA_model[["seg.q.tab"]] = SCNA_model_clonal_seg_Q_post( obs, b, delta, SCNA_model )
   seg.post_mix.w = SCNA_model[["seg.post_mix.w"]]
   
## Modeled seg info will be added to mode.tab
   SCNA_model[["seg.z.tab"]] = rowSums(seg.post_mix.w[, c("unif"), drop=FALSE]) 
   if( any(!is.finite(SCNA_model[["seg.z.tab"]]))) { stop() }

   SCNA_model[["seg.Wu.tab"]] = seg.post_mix.w[,"unif"] 
   SCNA_model[["seg.amp.tab"]] = seg.post_mix.w[,"amp"] 
   fac = (1-SCNA_model[["seg.z.tab"]])
   fac[ fac < 0] = 0  ## round error
   SCNA_model[["seg.qz.tab"]] = cbind( fac * SCNA_model[["seg.q.tab"]], SCNA_model[["seg.z.tab"]] )

## total CN if in allelic mode
   if( seg.obj[["copy_num_type"]] == "allelic" )
   {
      tot.obs = extract_total_copy_ratios_from_allelic_CAPSEG( seg.obj )  
      SCNA_model = apply_allelic_model_to_total_copy_ratios( tot.obs, b, delta, SCNA_model )
   }

  ## calculate allelic-balance, WGD
   if( seg.obj[["copy_num_type"]] == "allelic" )  
   {
     SCNA_model[["ab.tab"]] <- CalcAbDistr(obs, SCNA_model[["seg.qz.tab"]] )

     SCNA_model[["WGD"]] = ClassifySamplesWgdByProfile( seg.obj[["as.seg.dat"]], SCNA_model )
   }

   SCNA_model[["chr.arm.tab"]] <- CalcChrArmDistr(seg.obj, SCNA_model[["seg.q.tab"]], chr.arms.dat)

   return( SCNA_model )
}
