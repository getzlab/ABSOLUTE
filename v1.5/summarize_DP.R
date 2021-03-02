
summarize_mut_locations = function( DP_res, N_burn )
{
   N_iter = length(DP_res$c_fpost)
   N_mut = dim(DP_res$assign)[1]
   N_GRID = dim(DP_res$c_fpost[[2]])[2] 

#   loc = array( NA, dim=c(N_iter-N_burn+1, N_GRID, N_mut) )
   loc = array( 0, dim=c(N_GRID, N_mut) )

   for( i in N_burn:N_iter )
   {
      CL = unique( DP_res$assign[,i] )
      for( j in 1:length(CL) )
      {
         c_dat = which(DP_res$assign[,i] == CL[j] )
#         loc[ (i-N_burn+1), , c_dat ] = DP_res$c_fpost[[i]] [j,]
         loc[, c_dat ] = loc[, c_dat] + DP_res$c_fpost[[i]] [j,]
      }
   }

   res = t(loc / (N_iter-N_burn+1))

  ## should be N_mut X N_GRID
#   res = t( apply( loc, c(2,3), mean )  )

   return(res)
}



summarize_predictive_density = function( DP_res, N_burn, GRID, use_fixed=FALSE, fixed_dens=NA, fixed_scales=NA )
{
   N_iter = length(DP_res$c_fpost)
   pred = matrix( NA, nrow=N_iter-N_burn+1, ncol=length(GRID) )

   if( use_fixed )
   {
      fixed_mat = fixed_dens * fixed_scales
   }

   for( i in N_burn:N_iter )
   {
      if( use_fixed )
      {
         mat = rbind( DP_res$c_N[[i]] * DP_res$c_fpost[[i]], fixed_mat )
         pred[i-N_burn+1,] = colSums( mat ) / (sum(DP_res$c_N[[i]]) + sum(fixed_scales) )
      }
      else
      {
         pred[i-N_burn+1,] = colSums( DP_res$c_W[[i]] * DP_res$c_fpost[[i]] )
      }
   }

   res = apply( pred, 2, quantile, c(0.025, 0.5, 0.975) )
   res = rbind(res, colMeans(pred) )  

   return(res)
}


summarize_DP_num_components = function( DP_res, N_burn )
{
   DP_post = DP_res$DP_post
   DP_prior = DP_res$DP_prior
   k_0_map = DP_prior$k_0_map
   Pi_gamma = DP_prior$Pi_gamma

   N_mut = dim(DP_post$assign)[1]
   a = Pi_gamma$a  
   b = Pi_gamma$b

   N_iter = length(DP_post$c_fpost)
   est_K = rep( NA, N_iter-N_burn+1 )

   for( i in N_burn:N_iter )
   {
      est_K[(i-N_burn+1)] = DP_post$K[i]
   }

   return(est_K)
}


robust_summarize_DP_num_components = function( DP_res, N_burn, min_comp_size=1 )
{
   DP_post = DP_res$DP_post
   DP_prior = DP_res$DP_prior
   k_0_map = DP_prior$k_0_map
   Pi_gamma = DP_prior$Pi_gamma

   N_mut = dim(DP_post$assign)[1]
   a = Pi_gamma$a  
   b = Pi_gamma$b

   N_iter = length(DP_post$c_fpost)
   est_K = rep( NA, N_iter-N_burn+1 )

   for( i in N_burn:N_iter )
   {
#      est_K[(i-N_burn+1)] = DP_post$K[i]
      rr = rle( sort(DP_post[["assign"]][,i]))

      robust_K = sum(rr[["lengths"]] >= min_comp_size )
      est_K[(i-N_burn+1)] = robust_K
   }

   return(est_K)
}





## never used
cluster_DP_loc_results = function( DP_res, loc_res, GRID, N_burn )
{
  # identify modal # of DP components, and classify data accordingly
  
   est_k = summarize_DP_num_components( DP_res, N_burn )
   res=  rle( sort(est_k))
   use_k = res$values[  which.max(res$lengths) ]

   if( FALSE )
   {
    # version 0 - k-means on 1st moments of DP CCF
      E_loc = colSums( GRID * t(loc_res) )
      cres = kmeans( E_loc,  centers=use_k )
      assign = cres$cluster
      return(assign)
   }


   res = loc_clust_EM( loc_res, use_k, pi_clust=rep(5,use_k) )
   clust_p=res$clust_p   # N x K

   assign = apply( clust_p, 1, which.max )
   return(assign)
}



tree_cluster_DP = function( data, DP_res, N_burn, K, clust_CCF_prior, use_fixed=FALSE, fixed_dens=NA, fixed_scales=NA )
{
   DP_post = DP_res$DP_post

   N_iter = length(DP_post$c_fpost)
   N_mut = dim(DP_post$assign)[1]
   N_GRID = dim(DP_post$c_fpost[[2]])[2] 

   if( use_fixed )
   {
      N_fixed = nrow( fixed_dens )
      fixed_assign =  (N_mut+2):(N_mut+1+N_fixed)    ## assignment to these indices denotes fixed-cluster assignment
   }
   else { N_fixed = 0; fixed_assign = c() }

   coclust_mat = matrix( 0, nrow=N_mut+N_fixed, ncol=N_mut+N_fixed )

   for( i in N_burn:N_iter )
   {
      CL = unique( DP_post$assign[,i] )
      for( j in 1:length(CL) )
      {
         c_dat = sort( which(DP_post$assign[,i] == CL[j] ) )
         if( length(c_dat)== 1 ) { next } 

         for( k in 2:(length(c_dat)) ) 
         {
            for( l in 1:(k-1) ) 
            {
               coclust_mat[ c_dat[k], c_dat[l] ] = coclust_mat[ c_dat[k], c_dat[l] ]  + 1     
            }
         }
         if( CL[j] %in% fixed_assign )
         {
            fixed_ix = CL[j] - 1
            for( k in 1:(length(c_dat)) ) 
            {
               coclust_mat[ fixed_ix, c_dat[k] ] = coclust_mat[ fixed_ix, c_dat[k] ]  + 1     
            }
         }
      }
   }

   treeclust = hclust( as.dist(1/(coclust_mat+1)), method="complete" )
   assign = cutree( treeclust, k=K )
   CL = unique(assign)
   
   if( use_fixed )
   {
      fixed_mat = fixed_dens #* fixed_scales
      data = rbind(data, fixed_mat )
#      mat = rbind( DP_res$c_N[[i]] * DP_res$c_fpost[[i]], fixed_mat )
#      pred[i-N_burn+1,] = colSums( mat ) / (sum(DP_res$c_N[[i]]) + sum(fixed_scales) )
   }

   cluster_dens = matrix(NA, nrow=length(CL), ncol=dim(data)[2] )
   for( i in 1:length(CL))
   {
      c_ix = assign==CL[i]
      cluster_dens[i,] = get_cluster_stats( c_ix, log(t(data)), log(clust_CCF_prior) )
   }

   assign = assign[1:N_mut]
## TODO: sort clusters by CCF mode

   return( list("assign"=assign, "CCF_dens"=cluster_dens))
}


## initial version - does not handle fixed clusters
old_tree_cluster_DP = function( data, DP_res, N_burn, K, clust_CCF_prior )
{
   DP_post = DP_res$DP_post

   N_iter = length(DP_post$c_fpost)
   N_mut = dim(DP_post$assign)[1]
   N_GRID = dim(DP_post$c_fpost[[2]])[2] 

   coclust_mat = matrix( 0, nrow=N_mut, ncol=N_mut)

   for( i in N_burn:N_iter )
   {
      CL = unique( DP_post$assign[,i] )
      for( j in 1:length(CL) )
      {
         c_dat = sort( which(DP_post$assign[,i] == CL[j] ) )
         if( length(c_dat)== 1 ) { next } 

         for( k in 2:(length(c_dat)) ) 
         {
            for( l in 1:(k-1) ) 
            {
               coclust_mat[ c_dat[k], c_dat[l] ] = coclust_mat[ c_dat[k], c_dat[l] ]  + 1     
            }
         }
      }
   }

   treeclust = hclust( as.dist(1/(coclust_mat+1)), method="complete" )
   assign = cutree( treeclust, k=K )
   CL = unique(assign)
   
   cluster_dens = matrix(NA, nrow=length(CL), ncol=dim(data)[2] )
   for( i in 1:length(CL))
   {
      c_ix = assign==CL[i]
      cluster_dens[i,] = get_cluster_stats( c_ix, log(t(data)), log(clust_CCF_prior) )
   }

## TODO: sort clusters by CCF mode

   return( list("assign"=assign, "CCF_dens"=cluster_dens))
}





remove_empty_clusters = function( clust, N_DAT )
{
  dat_assign = clust[["assign"]][1:N_DAT]
  ps_assign = clust[["assign"]][-c(1:N_DAT)]

  ps_only = setdiff( unique(ps_assign), unique(dat_assign))

## delete clusters with no mutations assigned to them.  Can happen when pseudo-counts are used
   assign=dat_assign

   if( length(ps_only) > 0 ) 
   {
      d.ix = sort( ps_only, decreasing=TRUE )
      clust[["CCF_dens"]] = clust[["CCF_dens"]][-d.ix,,drop=FALSE]
   
      for( i in seq_along(d.ix) ) {
         assign[ assign > d.ix[i] ] = assign[ assign > d.ix[i] ] - 1
      }
   }

   clust[["assign"]] = assign

   return(clust)
}
