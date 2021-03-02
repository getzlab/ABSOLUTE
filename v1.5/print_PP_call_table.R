## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

PrintPpCallTable <- function(segobj.list, out.fn) {
  pp.call.tab <- data.frame()
  
  for (i in seq_along(segobj.list)) {
    segobj <- segobj.list[[i]]
    pp.vals <- GetPpTabVals(segobj, segobj[["mode.res"]][["call.status"]])
    pp.call.tab <- rbind(pp.call.tab, pp.vals )
  }

  write.table(pp.call.tab, file=out.fn, sep="\t", quote=FALSE, row.names=FALSE)
}

GetPpTabVals <- function(segobj, call.status) {
   val.names <- c("array", "sample", "call status", "purity", "ploidy",
                  "Genome doublings", "delta", "Coverage for 80% power",
                  "Cancer DNA fraction", "Subclonal genome fraction", "tau",
                  "E_CR")
   vals <- data.frame(matrix(NA, nrow=1, ncol=length(val.names)))
   colnames(vals) <- val.names

   vals[1, "array"] <- ifelse(is.null(segobj$array.name), segobj$sample.name, segobj[["array.name"]])
   vals[1, "sample"] <- segobj[["sample.name"]]

   mode.res <- segobj[["mode.res"]]
   
   if (call.status == "FAILED") {
     call.status <- mode.res[["mode.flag"]]
   }

   if (call.status %in% c("non-clonal", "non-aneuploid", "low purity")) {
     vals[1, "call status"] <- call.status
     return(vals)
   }

   a.s.d <- segobj[["as.seg.dat"]]
   e.cr <- sum(a.s.d[,"W"] * a.s.d[, "copy_num"])

   tau <- mode.res[["mode.tab"]][1, "tau"]
   gm <- mode.res[["mode.tab"]][1, "genome mass"]
   alpha <- mode.res[["mode.tab"]][1, "alpha"]
   wgd = mode.res[["mode.tab"]][1,"WGD"]
#   wgd <- ClassifySamplesWgdByProfile(segobj)
   subclonal.frac <- mode.res[["mode.tab"]][1, "frac.het"]
   delta <- alpha / (2*(1-alpha) + alpha*gm) 
   fDNA <- delta * gm

   SSNV_model = segobj[["mode.res"]][["mode_SSNV_models"]][[1]]
   pow80 <- NA
   if (delta > 0 & delta <= 1)
   {
      pow80 <- GetCovForPow(0.8, 1e-3/3, 0.5e-7, delta, SSNV_model)
   }
      
   vals[1, c(3:(ncol(vals)))] <- c(call.status, round(alpha,2), round(gm, 2),
                                   wgd, round(delta,2), pow80, round(fDNA, 2),
                                   round(subclonal.frac, 2), tau, e.cr)

   return(vals) 
}
