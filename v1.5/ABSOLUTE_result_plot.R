## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

## seg.dat is the result of ABSOLUTE

AbsoluteResultPlot <- function(sample.pdf.fn, seg.dat, called.mode.ix=NA,
                               verbose=FALSE) {
  pdf(sample.pdf.fn, 17.5, 18.5 )
    
  ## genome HSCR segplot
  d.mar <- par("mar")
  layout(mat=matrix(1:4, nrow=2, ncol=2, byrow=TRUE), widths=c(6, 2), heights=c(6, 6))

# par(mgp=)    # closer axis labels
# short axis ticks
  par(las=1)

#  if (!is.null(seg.dat[["allele.segs"]])) 
  if( TRUE )
  {
    allele.segs = get_hom_pairs_segtab(seg.dat)
    seg_colors = GetSegColsByAllelicBalance(seg.dat[["obs.scna"]], allele.segs  ) 
    PlotHscrAndSeghist(allele.segs, seg_colors, max_CR=2, plot.hist=TRUE )
    frame()
    frame()
  }
  par(mar=d.mar)



## Detailed plot for each mode found    
  ## PP mode plots
#  par(mfrow=c(4, 5))

# Single row per mode layout
#  layout( mat=matrix( 1:(4*7), nrow=4, ncol=7, byrow=TRUE), widths = c(1, 2.5, 1, 1, 1, 1, 1), heights=1 )

# Two rows of 6 plots each, for each mode
  PlotModes_layout()

# par(mgp=)    # closer axis labels
# short axis ticks

  PlotModes( seg.dat, n.print=nrow(seg.dat[["mode.res"]][["mode.tab"]]), 
             called.mode.ix=called.mode.ix, verbose=verbose)

  dev.off()
}
