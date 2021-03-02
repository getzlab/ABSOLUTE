## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

FilterSegs <- function(seg.info, min_probes=10, max_sd=100, verbose=FALSE) {  
  cn.mean <- mean(seg.info[, "copy_num"])
  cn.sd <- sd(seg.info[, "copy_num"])
  data.trunc <- cn.mean + as.numeric(max_sd) * cn.sd
  
  if (verbose) {
    print(paste(nrow(seg.info), " total segments", sep = ""))
  }
  
  rows <- (seg.info[, "n_probes"] >= min_probes) &
          (seg.info[, "copy_num"] <= data.trunc)
  
  seg.info <- seg.info[rows, ]
  seg.info <- seg.info[is.finite(seg.info[, "seg_sigma"]), ]
  
  if (verbose) {
    print(paste(nrow(seg.info), " segments remaining", sep = ""))
  }
  
  return(list(seg.info = seg.info))
}