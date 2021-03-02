## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

ApplyKaryotypeModel <- function(mode.res, model.id, train.obj, verbose=FALSE) {
    clust.res <- train.obj[[model.id]]$clust_res
    q <- dim(mode.res[["theta.q.tab"]])[2]
    m <- nrow(mode.res[["mode.tab"]])
    max.q <- dim(clust.res[["M_clust"]])[3]
    
    mode.res[["mode.clust.p"]] <- array(NA,
                                        dim = c(m, dim(clust.res[["M_clust"]])[1]))
    
    for (j in seq_len(m)) {
      mode.chrarm <- GetChrArmMat(mode.res, j, max.q)
      clust.ll <- DMultmix(mode.chrarm, clust.res, max.q, log.p=TRUE,
                           verbose=verbose)
      mode.res[["mode.tab"]][j, "Kar_LL"] <- clust.ll  
      mode.res[["mode.clust.p"]][j, ] <- MultmixClustP(mode.chrarm,
                                                       clust.res, max.q,
                                                       verbose=verbose)
    }
    
#    new.ll <- mode.res[["mode.tab"]][, "Kar_LL"] +
#              mode.res[["mode.tab"]][, "combined_LL"]
#    mode.res[["mode.tab"]][, "combined_LL"] <- new.ll
    
    return(mode.res)
}

GetChrArmMat <- function(mode.res, mode, max.q, n.arm= 39) 
{
  chr.arm.tab = mode.res[["mode_SCNA_models"]][[mode]][["chr.arm.tab"]]

  if (dim(chr.arm.tab)[2] == 2) {
    chrarm.mat <- rbind(chr.arm.tab[1, c(1:n.arm),],
                        chr.arm.tab[2, c(1:n.arm), ])
  } else {
    chrarm.mat <- chr.arm.tab[1, c(1:n.arm), ]
  }
  chrarm.mat <- chrarm.mat[, seq_len(max.q)]
  nas <- apply(is.na(chrarm.mat), 1, sum)
  rs <- rowSums(chrarm.mat)
  ix <- which(nas > 0 | rs != 1)
  chrarm.mat[ix, ] <- matrix(1 / max.q, nrow = 1, ncol = max.q, byrow = TRUE)
  
  return(chrarm.mat)
}

DMultmix <- function(theta, clust.res, max.q, log.p=FALSE, verbose=FALSE) {
  m.clust <- clust.res[["M_clust"]]
  clust.w <- clust.res[["clust_W"]]
  
  if (nrow(theta) == dim(m.clust)[2] / 2) {
    if (verbose) {
#      print("Switching to total copy Karyotype model...")
    }
    m.clust <- CollapseChrarmModelsToTcn(m.clust)
  } else {
    if (!(nrow(theta) == dim(m.clust)[2])) {
      stop("non-conforming arguments")
    }
  }
  
  ## disable outliers
  clust.w[length(clust.w)] <- 0
  clust.w <- clust.w / sum(clust.w)
    
  dat <- array(NA, dim = c(1, dim(theta)))
  dat[1, , ] <- theta
    
  mix.ll <- CalcDatLL(dat, m.clust, clust.w)
    
  if (log.p) {
    return(mix.ll)
  } else {
    return(exp(mix.ll))
  }
}

CollapseChrarmModelsToTcn <- function(m.clust) {
  k <- dim(m.clust)[1]  # clusters
  n.feat <- dim(m.clust)[2] / 2
  new.m.clust <- array(NA, dim = c(k, n.feat, dim(m.clust)[3]))
  
  for (k in 1:k) {
    for (i in 1:n.feat) {
      c_res = convolve(m.clust[k, i, ], rev(m.clust[k, i + n.feat, ]), type="open")
      c_res = c_res[c(1:(dim(m.clust)[3]))]
      c_res = c_res / sum(c_res)
      new.m.clust[k, i, ] = c_res
    }
  }
  
  return(new.m.clust)
}

MultmixClustP <- function(theta, clust.res, max.q, verbose=FALSE) {
  m.clust <- clust.res[["M_clust"]]
  clust.w <- clust.res[["clust_W"]]
  
  if (nrow(theta) == dim(m.clust)[2] / 2) {
    if (verbose) {
      print("Switching to total copy Karyotype model...")
    }
    m.clust <- CollapseChrarmModelsToTcn(m.clust)
  } else {
    if (!(nrow(theta) == dim(m.clust)[2])) {
      stop("non-conforming arguments")
    }
  }
  
  ## disable outliers
  clust.w[length(clust.w)] <- 0
  clust.w <- clust.w / sum(clust.w)
  
  dat <- array(NA, dim = c(1, dim(theta)))
  dat[1, , ] <- theta
  
  clust.p <- GetClustP(dat, m.clust, clust.w)
  
  ## sort clusters into plotting order
  if (!is.null(clust.res[["k_ix"]])) {
    clust.p <- clust.p[clust.res[["k_ix"]]]
  }
  
  return(clust.p)
}
