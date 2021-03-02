## These functions work for 1d or 2d (pairwise) CCF distributions (data)

run_DP = function( data, N_iter, Pi_k, clust_CCF_prior, use_fixed=FALSE, fixed_dens=NA, fixed_scales=NA )
{
   if(any(data < 0) ) { stop("Negative histogram bin!") }
   if(any( abs(rowSums(data) - 1) > 1e-10 )) { stop("Non-normalized histogram!") }

   N_dat = dim(data)[[1]]
   DP_prior = init_DP_prior(N_dat, Pi_k) 

   est = list()
   est$K = rep(NA, (N_iter+1))

   est$eta = rep(NA, (N_iter+1))
   est$gamma = rep(NA, (N_iter+1))
   est$gamma[1] = rgamma(1, DP_prior$Pi_gamma$a, DP_prior$Pi_gamma$b)

   est$assign = array(NA, dim=c(N_dat, (N_iter+1)) )
   est$assign[,1] = c(1:N_dat)
   est$c_fhat = list()
   est$c_W = list()
   est$c_N = list()
   est$c_fpost = list()

   use_gamma = 1
   for( i in 1:N_iter )
   {
      if( use_fixed)
      {
         cluster_res = DP_post_fixed( data, clust_CCF_prior, est$gamma[i], est$assign[,i], fixed_dens, fixed_scales )
      }
      else
      {
         cluster_res = DP_post( data, clust_CCF_prior, est$gamma[i], est$assign[,i])
      }

      est$assign[,(i+1)] = cluster_res$assign
      est$K[ (i+1) ] = length( unique(cluster_res$assign))

      est$eta[i+1] = rbeta( 1, est$gamma[i]+1, N_dat )   
      est$gamma[(i+1)] = sample_gamma_cond_N_k( N_dat, est$K[(i+1)], est$eta[(i+1)],  DP_prior$Pi_gamma ) 

      est$c_fhat[[(i+1)]] = cluster_res$c_fhat
      est$c_W[[(i+1)]] = cluster_res$c_W
      est$c_N[[(i+1)]] = cluster_res$c_N
      est$c_fpost[[(i+1)]] = cluster_res$c_fpost

      cat( est$K[i+1] )
      cat(" ")
   }

   return( list("DP_post"=est, "DP_prior"=DP_prior) )
} 


## DP API for histograms
cluster_prob = function( c_res, CCF_dist, n )
{
   Pr = sum( c_res * CCF_dist[n,] )
   return( Pr ) 
}

log_cluster_prob = function( log_c_res, log_CCF_dist, n )
{
   pr = log_c_res + log_CCF_dist[n,] 
   if( all(pr==-Inf) ) { return(-Inf) }
   else {
      log_Pr = LogAdd( pr )
   }
   return( log_Pr ) 
}


get_log_cluster_stats = function( c_ix, log_t_data, log_prior  )
{
   log_f_LL = rowSums( log_t_data[,c_ix,drop=FALSE]) + log_prior
   log_f_post = log_f_LL - LogAdd(log_f_LL)   

   return(log_f_post)
}


old_get_cluster_stats = function( c_ix, data  )
{
   log_f_LL = colSums( log(data[c_ix,,drop=FALSE]) )
   log_f_post = log_f_LL - LogAdd(log_f_LL)   
   f_post = exp(log_f_post) 

   return(f_post)
}

get_cluster_stats = function( c_ix, log_t_data, log_prior  )
{
   log_f_LL = rowSums( log_t_data[,c_ix,drop=FALSE]) + log_prior
   log_f_post = log_f_LL - LogAdd(log_f_LL)   
   f_post = exp(log_f_post)

   return(f_post)
}

