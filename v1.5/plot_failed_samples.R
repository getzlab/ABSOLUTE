## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

PlotFailedSamples <- function(segobj.list, pdf.fn) {
  if (length(segobj.list) == 0) {
    return()
  }
  
  Q <- 15
  x.max <- 2.5

  ## loop over samples to get range of seg std errs and MODE_FLAG values
  stderrs <- c()
  mode.flags <- rep(NA, length(segobj.list))
  
  for (i in 1:length(segobj.list)) {
    obs <- ExtractSampleObs(segobj.list[[i]])

    seg_sigma_num = 0.5  ## TODO - get rid of this - not used in downstream model - but crashes filtering code if missing
    stderrs = c(stderrs, seg_sigma_num / sqrt(as.numeric(obs[["n_probes"]])))


    mode.flags[i] <- segobj.list[[i]][["mode.res"]][["mode.flag"]]
  }
  
  ## cap at 95th percentile
  p95 <- sort(stderrs)[floor(length(stderrs) * 0.95)]
  stderrs[stderrs > p95] <- p95
  
  ## sort samples by MODE_FLAG
  ix <- sort(mode.flags, index.return=TRUE)[["ix"]]
  segobj.list <- segobj.list[ix]

  pdf(pdf.fn, 14, 15)
  par(mfrow=c(5, 4))
  pal <- colorRampPalette(c("blue", "purple", "lightgray"))
  colpal <- pal(1000)
  
  for (i in seq_along(segobj.list)) {
    obs <- ExtractSampleObs(segobj.list[[i]])
    d <- obs[["d"]]
    ix <- (d >= 0) & (d < x.max)
    d.stderr <- obs[["d.stderr"]]
    d.stderr[d.stderr > p95] <- p95
    
    PlotSeglenHist(d[ix], obs[["W"]][ix], log(d.stderr[ix]),
                   color.range=range(log(stderrs)), use.pal=colpal,
                   bin.w=0.025, x.max=x.max)
    title(sub=segobj.list[[i]][["array"]])
    mtext(segobj.list[[i]][["sample.name"]], line=1, cex=0.75, adj=0)
    
    hscn_params <- obs[["error.model"]]
    if (!is.null(hscn_params$sigma.eta)) {
    mtext(substitute(paste(sigma[eta] == x, ",", sigma[epsilon]==y, sep=""),
                     list(x=round(hscn_params[["sigma.eta"]], 3),
                          y=round(hscn_params[["sigma.nu"]], 3))),
          line=0, adj=1, cex=0.75)
    }
    mtext(segobj.list[[i]][["mode.res"]][["mode.flag"]], line=1, adj=1, cex=0.75)
  }
  
  dev.off()
}

PrintFailedTable <- function(segobj.list, failed.tab.fn) {
  tab <- matrix(NA, nrow=length(segobj.list), ncol=3)
  colnames(tab) <- c("array", "sample", "mode.flag")
  
  for (i in seq_along(segobj.list)) {
    tab[i, "array"] <- segobj.list[[i]][["array.name"]]
    tab[i, "sample"] <- segobj.list[[i]][["sample.name"]]
    tab[i, "mode.flag"] <- segobj.list[[i]][["mode.res"]][["mode.flag"]]
  }
  
  write.table(tab, file=failed.tab.fn, row.names=FALSE, sep="\t", quote=FALSE)
}

