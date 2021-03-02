## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CalcDatLL <- function(dat, m.clust, clust.w) {
  k <- dim(m.clust)[1]
  clust.ll <- GetDatClustLL(dat, m.clust)
  
  clust.w <- matrix(clust.w, ncol = (k + 1),
                    nrow = nrow(clust.ll), byrow = TRUE)
  clust.ll <- clust.ll + log(clust.w)    
  return(LogAdd(clust.ll))
}

GetDatClustLL <- function(dat, m.clust) {
  ## dat: samples X features X Q 
  N <- dim(dat)[1]
  Q <- dim(dat)[3]
  k <- dim(m.clust)[1]
  J <- dim(m.clust)[2]
  
  feat.clust.ll <- array(NA, dim = c(N, J, k + 1))
  
  ## clust_LL:  samples X clusters   
  ## p^x
  for (k in 1:k) {
    for (j in 1:J) {
      log.p <- matrix(log(m.clust[k, j, ]), nrow = N, ncol = Q, byrow = TRUE)
      feat.clust.ll[, j, k] <- rowSums(log.p * dat[, j, ])
    }
  }
  
  ## outlier cluster
  for (j in 1:J) {
    log.p <- matrix(log(rep(1/Q, Q)), nrow = N, ncol = Q, byrow = TRUE)
    feat.clust.ll[, j, k + 1] <- rowSums(log.p * dat[, j, ])
  }
  
  ## sum over all features
  clust.ll <- apply(feat.clust.ll, c(1, 3), sum)
  
  return(clust.ll)
}

GetClustP <- function(dat, m.clust, clust.w) {
  k <- dim(m.clust)[1]
  clust.ll <- GetDatClustLL(dat, m.clust)
  
  clust.w <- matrix(clust.w, ncol = (k + 1), nrow = nrow(clust.ll), byrow = TRUE)
  clust.ll <- clust.ll + log(clust.w)
  
  clust.p <- exp(clust.ll - LogAdd(clust.ll))
  
  return(clust.p)
}
