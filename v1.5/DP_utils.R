get_log_stirling_coefs = function( W )
{
 ## if W = c(1:N), then this function returns the (unisgned 1st kind of) Stirling numbers for
 ## N, k=c(1:N)
 ## starts to give incorrect results at N = 19 if fft is used

## x and y logged
   log_conv = function( x, y )
   {
       ## y is len 2
      x = c(-Inf,x,-Inf)

      N = length(x) - 2
      res = rep(NA,N+1)

      for( k in 1:(N+1) )   
      {
         res[k] = LogAdd(  c(x[k] + y[1],  x[k+1] + y[2]) ) 
      }
      return(res)
   }


   N = length(W)
   nW = W 

   cres = log(1)
   for( i in 1:(N-1) )
   {
      cres = log_conv( cres, log(rev( c(1, nW[i])) ) ) 
   }

   return( rev(cres) )
}




get_stirling_coefs = function( W )
{
 ## if W = c(1:N), then this function returns the (unisgned 1st kind of) Stirling numbers for
 ## N, k=c(1:N)
 ## starts to give incorrect results at N = 19 if fft is used

 ##  this implemtation overflows at N=141.    Use the log impementation instead

   my_conv = function( x, y )
   {
       ## y is len 2
      x = c(0,x,0)

      N = length(x) - 2
      res = rep(0,N+1)

      for( k in 1:(N+1) )   
      {
#         res[k] = x[k-1 +1] * y[1] + x[k-1 +2] * y[2] 

         res[k] = x[k] * y[1] + x[k+1] * y[2] 
      }
      return(res)
   }


   N = length(W)
#   nW = N * W 
   nW = W 

#   cres = c(0,nW[1])
   cres = 1
   for( i in 1:(N-1) )
   {
#      cres = convolve( cres, rev( c(1, nW[i]) ), type="o" ) 
      cres = my_conv( cres, rev( c(1, nW[i]) ) ) 
   }

   return( rev(cres) )
}

# map 1st moments of (k | N, gamma), over a grid on gamma
get_k_0_map = function( N, gamma_GRID, log_stirling_coef )
{
#   gamma_GRID = seq(0.01,10, length.out = 1000) 
   N_gamma = length(gamma_GRID)

   k_0_map = matrix( 0, ncol=2, nrow=N_gamma )
   colnames( k_0_map ) = c("k_0", "gamma" )
   k_0_map[,2] = gamma_GRID 

   DP_prob = array( 0, dim=c(N_gamma, N) )
   for(i in 1:N_gamma )
   {
      DP_prob[i,] = DP_prob_k_cond_gamma_N( N, gamma_GRID[i], log_stirling_coef )
      k_0_map[i,1] = sum( DP_prob[i,] * c(1:N) )
   }

   return( k_0_map )
}

get_gamma_prior_from_k_prior = function( N, k_0_map, k_prior )
{
   LL = function( Par, N, int_k_prior )
   {
      if( any(Par < 0) ) { return(Inf) }

      mu = Par[1]
      sigma = Par[2]
      B = mu/sigma^2
      A = B*mu

      k_prob = dgamma( k_0_map[,2],  A, B )
      k_prob = k_prob * c( 1, 1/diff(c( k_0_map[,1]) ) )
      k_prob = k_prob / sum(k_prob)

      ix = k_prob > 0 & int_k_prior>0
      t1 = sum( (log(k_prob)*k_prob)[ix] )  
      t2 = sum( (log(int_k_prior)*k_prob)[ix] )
      divergence =  (t1 - t2) / length(int_k_prior) 

      return( divergence )
   }    

   sigma_grid = seq(1, 25, by=5 )
   mu_grid = seq(1, 25, by=5 )
   obj = array(NA, dim=c(length(sigma_grid), length(mu_grid)) )
   mode_vals = array( NA, dim=c( length(sigma_grid)*length(mu_grid), 3 ) )

   int_k_prior = approx( x=c(1:N), y=k_prior, xout=k_0_map[,1])$y 
#   int_k_prior = int_k_prior * c(1, 1/diff(c(k_0_map[,1])) )
   int_k_prior = int_k_prior / sum(int_k_prior)

   for( i in 1:length(sigma_grid) )
   {
      for( j in 1:length(mu_grid))
      { 
         par0=c(sigma_grid[i], mu_grid[j])
         res = optim( par0, fn=LL, method="Nelder-Mead", N=N, int_k_prior=int_k_prior )
         val = res$par
         obj[i,j] = res$value
         mode_vals[ (i-1)*length(mu_grid) + j , ] = c(res$par, res$value)
        cat('.')
      }
   }
   cat("\n")

   ix = which.min( mode_vals[,3])
   mu = mode_vals[ix,1]
   sigma = mode_vals[ix,2] 
   KL_divergence = mode_vals[ix,3]

   B = mu/sigma^2
   A = B*mu
   val = list("a"=A, "b"=B, "KL_divergence"=KL_divergence )

   return(val)
}


## Pi_k parameterizes a negative binomial prior distribution over the number of DP components (k).
## This function returns parameters of a prior gamma distribution specifying the best fit to this distribution, taking into account the number of data elements (N).
init_DP_prior = function(N, Pi_k)
{
#   p = mu / (r*(1+mu))
   k_prior = dnbinom( c(1:N), size=Pi_k$r, mu=Pi_k$mu )
   k_prior = k_prior / sum(k_prior) 

   print( paste( "Initializing prior over DP k for ", N, " items", sep="") )

   log_stirling_coef = get_log_stirling_coefs(c(1:N)) 
   GMAX= 5
   GRID = seq(1e-25, GMAX, length=1000)
   k_0_map = get_k_0_map( N, GRID, log_stirling_coef )

   Pi_gamma = get_gamma_prior_from_k_prior( N, k_0_map, k_prior )

   return(list("log_stirling_coef"=log_stirling_coef, "Pi_gamma"=Pi_gamma, "Pi_k"=Pi_k, "k_prior"=k_prior, "k_0_map"=k_0_map) )
}


test_prior = function(N, lambda, mu)
{
   par( mfcol=c(2,2) )
   par( las=1 )

   DP_prior = init_DP_prior(N, Pi_k=list("r"=lambda, "mu"=mu ) )
   k_0_map = DP_prior$k_0_map
   k_prior = DP_prior$k_prior
   XMAX = qnbinom( p=0.99,  size=lambda, mu=mu )

   print(  pgamma( max(k_0_map[,2]),  DP_prior$Pi_gamma$a, DP_prior$Pi_gamma$b, lower.tail=TRUE)  )
   DG = dgamma( k_0_map[,2],  DP_prior$Pi_gamma$a, DP_prior$Pi_gamma$b ) 
   k_prob = DG
   k_prob = k_prob/sum(k_prob) * c( 1, 1/diff(c(k_0_map[,1])) )

   YLIM = max( c(k_prior, k_prob) )
   plot( 0, type='n', bty="n",  xlab="Number of DP components (K)", ylab="Density", main="", ylim=c(0,YLIM), xlim=c(0, XMAX) )

   min_K = k_0_map[1,1]
   int_k_prior = approx( x=c(1:N), y=k_prior, xout=k_0_map[,1])$y 
#   int_k_prior = int_k_prior/sum(int_k_prior) * c(1, 1/diff(c(k_0_map[,1])) ) 

   lines( k_0_map[,1], k_prob, col=1, lty=2 )
#   lines( k_0_map[,1], int_k_prior, col=3 )
   lines( c(1:N), DP_prior$k_prior, col=1 )

   print( sum(int_k_prior) )
   print( sum(k_prob) )
   print(DP_prior$Pi_gamma) 

   plot( 0, type='n', bty="n",  xlab="gamma", ylab="Density", main="", ylim=c(0,max(DG) ), xlim=c(0.01, max(k_0_map[,2])) )
   lines( k_0_map[,2], DG )
}

DP_prob_k_cond_gamma_N = function( N, gamma, log_stirling_coef )
{
 ## Escobar and West 1995 for DP loglik eqn 10.  
  ## note we change notation for gamma to alpha to follow convention
   alpha=gamma
   loglik = rep(NA, N)
  
   for( k in 1:N )
   {
      loglik[k] = log_stirling_coef[k] + lgamma(N-1) + k*log(alpha) + lgamma(alpha) - lgamma(alpha + N)
   }

   Pr = exp( loglik - LogAdd(loglik) )

   return(Pr)
}


sample_gamma_cond_N_k = function( N, k, eta, Pi_gamma )
{
  ###  Escobar and West 1995
   a = Pi_gamma$a
   b = Pi_gamma$b

   m1 = rgamma( 1, a + k, b - log(eta) ) 
   m2 = rgamma( 1, a + k - 1, b - log(eta) )

   D = N * (b - log(eta))
   w = (a+k-1) / (D+a+k)
   new_gamma = w * m1 + (1-w) * m2

   return(new_gamma)
}



get_perm_stirling_coefs = function( W )
{
   N = length(W)
   sum_coefs = rep(0, N)
   last_coefs = rep(1,N)
#   epsilon = 1e-6
   epsilon = 5e-4
   C = Inf
   i = 1

   while( C > epsilon )
   {
      conf_W = sample( W, length(W) )
      new_coefs = get_stirling_coefs( cumsum(conf_W) )
      sum_coefs = sum_coefs + new_coefs
      coefs = sum_coefs / i

      i = i + 1
      C = abs( 1 - sum( (coefs / last_coefs), na.rm=TRUE) / N )
  
      last_coefs = coefs

#      if( !is.finite(C) ) { C=epsilon+1; next }

cat("."); if( i %% 25 == 0 ) { print(C) }

   }
cat("\n")

   return( coefs )
}



alt_get_perm_stirling_coefs = function( W )
{
   N = length(W)
   sum_coefs = rep(0, N)
   last_coefs = rep(1,N)

   SET = 100
   coefs = matrix(0, nrow=SET, ncol=N )

   epsilon = 1e-2
   C = Inf
   i = 1; j=1

   while( C > epsilon )
   {
      conf_W = sample( W, length(W) )
      coefs[i,] = get_stirling_coefs( cumsum(conf_W) )
      i = i + 1

      if( i %% SET == 0 )
      {
         if( j > 1 )
         {
#            mm = apply( as.matrix(means[1:(j-1),]), 2, mean )
#            m_new = apply( means, 2, mean )
            
            lcoefs = log(coefs)
 
            mm = apply( as.matrix(lcoefs[1:(i-SET),]), 2, mean )
            m_new = apply( lcoefs[1:(i-1),], 2, mean )
            SDs = apply( lcoefs[1:(i-1),], 2, sd )

            vals = abs( m_new - mm ) / SDs 
            C = sum( vals[1:(N-1)] )
         }

         print(C)

         coefs = rbind( coefs, matrix( 0, nrow=SET, ncol=N) )
         j = j + 1
      }
      cat(".")

   }
cat("\n")

#   means = exp(  apply( log(coefs[(1:(i-1)),]), 2, mean )   )
   means =  apply( coefs[(1:(i-1)),], 2, mean )   
   return( means )
}


DP_loglik = function( N, k, gamma )
{

  ## Carter and Getz 2008 with seglen correction.
#   loglik = log(coefs[K]) + K*log(alpha) - log( eval_M_polynom(alpha, coefs) )

   return(loglik)
}



eval_M_polynom = function( alpha, coefs )
{
   pows = c(1:length(coefs))
   res = sum( coefs * alpha^pows  )

   return(res)
}
