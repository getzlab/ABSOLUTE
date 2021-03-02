## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GetCovForPow <- function(required.power, eps, fdr, delta, SSNV_model) 
{
  max <- 1000
  cov <- 2
  
  if (delta <= 0 | delta > 1) {
    return(NA)
  }
  
  while (TRUE) {
    pow <- PowerCalc(cov, eps, fdr, delta, SSNV_model)
    if (pow >= required.power) {
      break
    }
    if (cov > max) {
      cov <- NA
      break
    }
    
    cov <- cov + 1
  }
  return(cov)
}


# Carter et al 2012
orig_PowerCalc <- function(n, eps, fdr, delta) {
  p <- dbinom(c(0:n), size = n, prob = eps)
  pv <- 1 - cumsum(c(0, p))
  ks <- min(which(pv <= fdr))
  
  cval <- (fdr - pv[ks]) / (pv[ks - 1] - pv[ks])
  
  p1 <- dbinom(c(0:n), size = n, prob = delta)
  pow <- 1 - sum(p1[c(1:(ks - 1))]) + cval * p1[ks - 1]
  
  return(pow)
}


# new version takes SSNV skew into account
PowerCalc <- function(n, eps, fdr, delta, SSNV_model)
{
   f_skew = SSNV_model[["SSNV_skew"]]
   rho = SSNV_model[["rho"]]

 #  A = f_skew * eps * rho
 #  B = (1-f_skew * eps) * rho
 #  p <- d_beta_binom( c(0:n), A, B, n )
## probability of observing sequencing errors is not affected by allelic skew
   p <- dbinom(c(0:n), size = n, prob = eps)

   pv <- 1 - cumsum(c(0, p))
   ks <- min(which(pv <= fdr))
  
   cval <- (fdr - pv[ks]) / (pv[ks - 1] - pv[ks])

   A = f_skew * delta * rho
   B = (1-f_skew * delta) * rho
  
   p1 <- d_beta_binom(c(0:n), A, B, n )
   pow <- 1 - sum(p1[c(1:(ks - 1))]) + cval * p1[ks - 1]
  
   return(pow)
}

PowerCalc_for_single_read = function( n, delta, SSNV_model )
{
   f_skew = SSNV_model[["SSNV_skew"]]
   rho = SSNV_model[["rho"]]

   A = f_skew * delta * rho
   B = (1-f_skew * delta) * rho
  
   p1 <- d_beta_binom(0, A, B, n ) ## prob of observing 0 reads
   pow <- 1 - p1
#   if( pow < 0 ) { pow = 0 }  ## round error
  
   return(pow)

}

mode_SSNV_pow_calc = function( SSNV_model, mut.cn.dat, mut.modeled.cn, alpha, single_read=FALSE )
{
   N = nrow(mut.cn.dat)
   pow = rep(NA, N)

## appropriate for mutect without forced-calling
## TODO: expose these at top level
   eps = 1e-3/3
   fdr = 0.5e-7

   cov = rowSums(mut.cn.dat[,c("alt", "ref")])
   qt = mut.modeled.cn[,"q_hat"]
## allelic fraction of CCF=1 mutation at mult 1
   delta = alpha / (2*(1-alpha) + alpha * qt )

   if( single_read == FALSE )
   {
      for( i in 1:N )
      {
         pow[i] = PowerCalc( cov[i], eps, fdr, delta[i], SSNV_model )
      }
   }
   else
   {
      for( i in 1:N )
      {
         pow[i] = PowerCalc_for_single_read( cov[i], delta[i], SSNV_model )
      }
   }

   return(pow)
}

