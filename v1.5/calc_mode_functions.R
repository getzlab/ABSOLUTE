## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CalcModeLogCurv <- function(mode, mode.hess, verbose=FALSE) {
   hess.mat <- mode.hess
   b <- mode[1]

   curvature <- abs(det(hess.mat / (2 * pi)))

   mode.curv <- log((curvature)^(-1/2) ) - log(2)

   if ((!is.finite(mode.curv)) && verbose) {
     print("WARNING: NON-FINITE log_evidence")
   }
   
   return(mode.curv)
}

CalcAbDistr <- function(obs, seg.qz) {
   uniq.seg.ix <- unique(obs[["seg.ix"]])
   n.seg <- length(uniq.seg.ix) 
   seg.ab.tab <- array(NA, dim=c( n.seg, ncol(seg.qz)))

   found.pairs <- rep(FALSE, n.seg)
   use.w <- rep(NA, n.seg)
   for (i in 1:n.seg) {
     if (sum(obs[["seg.ix"]] == uniq.seg.ix[i]) != 2 ) {
       next
     }
     
     found.pairs[i] <- TRUE
     pair.ix <- which(obs[["seg.ix"]] == uniq.seg.ix[i]) 
     
     use.w[i] <- obs[["W"]][pair.ix[1]] 
     seg.ab.tab[i, ] <- apply(seg.qz[ pair.ix, ], 2, prod)
   }
   
   seg.ab.tab <- seg.ab.tab[found.pairs, ]
   use.w <- use.w[found.pairs]
   use.w <- use.w / sum(use.w)
   ab.tab <- colSums(use.w * seg.ab.tab)
   return(ab.tab)
}

CalcFitEntropy <- function(obs, seg.qz) {
  H <- function(p) {
    v <- -p * log(p)
    v[p == 0] <- 0
    return(sum(v))
  }

  seg.h <- apply(seg.qz, 1, H)
  res <- sum(obs[["W"]] * seg.h )
  
  return(res)
}


