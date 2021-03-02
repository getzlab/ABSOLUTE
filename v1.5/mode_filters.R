## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GenomeHetFilter <- function(obs, mode.res, max.non.clonal, max.neg.genome,
                            Q, verbose=FALSE) {
  ## calculate provisional seg_Z_tab and filter out modes that imply > 50% genome het.
  ## and filter out modes with > 2.5% het genome < 0
  mode.tab <- mode.res[["mode.tab"]]
  ## init both to all zeros
  frac.het <- frac.neg.het <- rep(0, nrow(mode.tab))

  for (i in seq_len(nrow(mode.tab))) {
    delta <- mode.tab[i, "delta"]
    b <- mode.tab[i, "b"]
    
    obs[["error.model"]][["fit.at"]] <- mode.tab[1, "AT"]
    comb <-  GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
#    seg.z <- mode.res[["seg.qz.tab"]][i, , Q+1]
    seg.z <- mode.res[["mode_SCNA_models"]][[i]][["seg.qz.tab"]][,Q+1]

    frac.het[i] <- sum(seg.z * obs[["W"]])
#    frac.neg.het[i] <- sum((obs[["W"]] * seg.z)[obs[["d.tx"]] < comb[1]])
    frac.neg.het[i] <- sum((obs[["W"]])[ seg.z > 0.9 & obs[["d.tx"]] < comb[1]])
  }

  if (max.non.clonal > 0) {
    nc.ix <- (frac.het > max.non.clonal)
  
    if (verbose) {
      print(paste("removing ", sum(nc.ix), " / ", length(nc.ix),
                  " modes with >", max.non.clonal*100, "% genome non-clonal.", sep=""))
    }
  } else {
    nc.ix <- rep(FALSE, length(frac.het))
  }

if( FALSE )
{
  if (max.neg.genome > 0) {
    neg.mode.ix <- (frac.neg.het > max.neg.genome) & (!nc.ix)

  if( is.na(neg.mode.ix)){ stop() }

    if (verbose) {
      print(paste("removing ", sum(neg.mode.ix), " / ", length(neg.mode.ix),
                  " modes with >", (max.neg.genome * 100) ,
                  "% genome non-clonal < 0 copies.", sep=""))
    }
  } else {
    neg.mode.ix <- rep(FALSE, length(frac.het))
  }
}

  if( any(is.na(nc.ix)) ) { stop() }

  ## return the 'bad' indices
  return(nc.ix)

}





NegGenomeFilter = function( obs, mode.tab, max.neg.genome, Q, verbose=FALSE) 
{
  frac.neg.het <- rep(0, nrow(mode.tab))
  eps = 0.1 ## this much below comb(0)

  for (i in seq_len(nrow(mode.tab))) {
    delta <- mode.tab[i, "delta"]
    b <- mode.tab[i, "b"]
    
    obs[["error.model"]][["fit.at"]] <- mode.tab[1, "AT"]
## TODO: add X-chr support
    comb_A <-  GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
#    frac.neg.het[i] <- sum((obs[["W"]] * seg.z)[obs[["d.tx"]] < comb[1]])

    frac.neg.het[i] <- sum( obs[["W"]][ obs[["d.tx"]] < comb_A[1] - eps ] )
  }


  if (max.neg.genome > 0) {
    neg.mode.ix <- (frac.neg.het > max.neg.genome)

  if( is.na(neg.mode.ix)){ stop() }

    if (verbose) {
      print(paste("removing ", sum(neg.mode.ix), " / ", length(neg.mode.ix),
                  " modes with >", (max.neg.genome * 100) ,
                  "% genome non-clonal < 0 copies.", sep=""))
    }
  } else {
    neg.mode.ix <- rep(FALSE, length(frac.neg.het))
  }

  return( neg.mode.ix )
}
