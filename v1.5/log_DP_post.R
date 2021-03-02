DP_post = function( data, prior, gamma, assign_init )
{
   log_t_data = log(t(data))  ## for speed optimization
   log_data = log(data)  
   log_prior = log(prior)

   N = nrow(data)
   assign = assign_init
   clusters = array(NA, dim=N+1)
   cluster_stats = list(N+1)

   for( j in 1:(N+1) )
   {
      cluster_stats[[j]] = 0
      clusters[j] = sum(assign == j)
   }

   c_loglik = array(-Inf, dim=N+1)
   dirty = array(TRUE, dim=N+1)
   dirty[ clusters == 0 ] = FALSE

   for( n in 1:N )
   {
      c_loglik = rep(-Inf, length(c_loglik) )
      skip_count = 1

      NZ_clusters = which( clusters != 0 )
      for( j in 1:length(NZ_clusters) )
      {
         c_ix = NZ_clusters[j]
         c_dat = (assign == c_ix)

         if( dirty[c_ix] == TRUE )
         {
            res = get_log_cluster_stats( c_dat, log_t_data, log_prior )
            cluster_stats[[c_ix]] = res
            dirty[c_ix] = FALSE
         }

       ## if the current point is the only thing in the cluster...
        # This seems to work empirically (as well as theoretically)
         if( sum(c_dat) == 1 & c_dat[n] == TRUE )
         {
            skip_count = skip_count + 1   ## at most = 2
            next
         }

         const = sum(c_dat) / (N - 1 + gamma)

         if( c_dat[n] == FALSE )
         {
            c_res = cluster_stats[[c_ix]]
         }
         else
         {
            ## remove the current point from the cluster
            ## and re-compute cluster stats

            c_dat[n] = FALSE
            c_res = get_log_cluster_stats( c_dat, log_t_data, log_prior )
         }

         log_c_pr = log_cluster_prob(  c_res, log_data, n  )
         c_loglik[c_ix] = log(const) + log_c_pr
      }

    ## prob to open a new cluster
      const = (gamma) / (N - 1 + gamma) 

#      temp_c_dat = rep(FALSE, length(c_dat) )
#      temp_c_dat[n] = TRUE
#      new_c_res = get_log_cluster_stats( temp_c_dat, log_t_data, log_prior )
#      log_c_pr = log_cluster_prob( new_c_res, log_data, n )

      new_c_res = log_prior
      log_c_pr = log_cluster_prob( new_c_res, log_data, n )

      h_loglik = log(skip_count) + log(const) + log_c_pr

      new_c_ix = min( which(clusters==0) )
      c_loglik[new_c_ix] = h_loglik

      c_loglik[ is.nan(c_loglik)] = -Inf
      c_lik = exp(c_loglik - LogAdd(c_loglik))

      if( sum( c_lik ) > 0 )
      {
         res = which( rmultinom( 1, 1, c_lik ) == 1 )
      }
      else
      {
         print("DP: all c_lik==0!!!")
         res = new_c_ix
      }

     ## if a point has been assigned to a different cluster..
      if( assign[n] != res )
      {
         clusters[ assign[n] ] = clusters[ assign[n] ] - 1
         if( clusters[ assign[n] ] > 0 )
         {
            dirty[assign[n] ] = TRUE
         } 
         else
         {
            cluster_stats[[ assign[n] ]] = 0
         }
         dirty[res] = TRUE
         clusters[res] = clusters[res] + 1
         assign[n] = res
      }
   }

   NZ_clusters = unique(assign)

   c_fhat = rep(0, length(NZ_clusters) )
   c_W = rep(0, length(NZ_clusters) )
   c_N = rep(0, length(NZ_clusters) )
   c_fpost = matrix(NA, nrow=length(NZ_clusters), ncol=ncol(data) )

   for( j in 1:length(NZ_clusters) )
   {
      c_ix = NZ_clusters[j]
      c_dat = (assign == c_ix)

      if( dirty[c_ix] == TRUE )
      {
         res = get_log_cluster_stats( c_dat, log_t_data, log_prior )
         cluster_stats[[c_ix]] = res
         dirty[c_ix] = FALSE
      }

      c_res = exp(cluster_stats[[c_ix]])

      c_fhat[j] = which.max(c_res) 
      c_W[j] = sum(c_dat) / length(c_dat) 
      c_N[j] = sum(c_dat)
      c_fpost[j,] = c_res
   }

#   DP_loglik = calc_DP_loglik( N, length(NZ_clusters), gamma, stirling_coef )
   DP_loglik = NA

   return( list("c_fhat"=c_fhat, "c_W"=c_W, "c_N"=c_N, "c_fpost"=c_fpost, "assign"=assign, "loglik"=DP_loglik, "K"=length(NZ_clusters) ) )
}

