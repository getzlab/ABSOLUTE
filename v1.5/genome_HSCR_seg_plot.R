## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

get_genome_coords = function( chrpos )
{ 
  chr.lens <- GetChrLens(x=FALSE)
  chr.lens <- as.numeric(chr.lens)
  chr.w <- chr.lens / sum(chr.lens)
  chr.offsets <- c(0, cumsum(chr.w[c(1:(length(chr.w) - 1))]), 1)

  genome.coords = rep(NA, nrow(chrpos) )
  for( i in 1:nrow(chrpos) )
  {
    chr <- as.integer( chrpos[i,1] )
    genome.coords[i] = chr.offsets[chr] + as.numeric(chrpos[i,2]) / chr.lens[chr] * chr.w[chr]
  }

  return( genome.coords )
}

GenomeHscrSegPlot <- function(allele.segs, seg_colors, y.lab, y.min, y.max, plot.model.fit=FALSE, plot.comb.fit=FALSE, mode.info=NA, comb=NA, comb.color=NA, plot.seg.sem=FALSE, plot.total.CN=FALSE, total.seg.dat=NA, tot.seg.colors=NA, log2CR = FALSE, label.chrs=TRUE ) 
{
  chr.lens <- GetChrLens(x=FALSE)
  chr.lens <- as.numeric(chr.lens)
  chr.w <- chr.lens / sum(chr.lens)
  
  plot(0, type = "n", bty = "n", xlim = c(0, 1), ylim = c(y.min, y.max),
       xlab = "", ylab = y.lab, main = "", xaxt = "n")

#  mtext( "Chromosome", side=1, line=1, adj=0 )
  
  if( label.chrs )
  {
     ww <- as.vector(rbind(chr.w, chr.w)) / 2
     chr.mids <- cumsum(ww)[(c(1:length(ww)) - 1) %% 2 == 0]
    
     lab.vals <- (c(1:length(chr.w)))
     odd.ix <- lab.vals %% 2 == 1
    
     mtext(text = lab.vals[odd.ix], side = 1, line = -0.45, at = chr.mids[odd.ix], 
           las = 1, cex = par("cex.axis") * par("cex") * 0.9)
     mtext(text = lab.vals[!odd.ix], side = 1, line = 0, at = chr.mids[!odd.ix],
           las = 1, cex = par("cex.axis") * par("cex") * 0.9)
   }
    
  chr.offsets <- c(0, cumsum(chr.w[c(1:(length(chr.w) - 1))]), 1)
  cent.pos <- (GetCentromerePos())/sum(chr.lens) +
    chr.offsets[c(1:(length(chr.offsets) - 1))]
    
  for (i in 1:(length(chr.offsets) - 1)) {
    use.col <- ifelse(i%%2 == 1, "grey90", "white")
    rect(xleft = chr.offsets[i], ybottom = y.min, xright = chr.offsets[i + 1], 
         ytop = y.max, col = use.col, border = NA)
    lines(y = c(y.min, y.max), x = rep(cent.pos[i], 2), lty = 3, lwd = 0.5)
  }

  AS.seg.ix = allele.segs[, c("seg.ix.1", "seg.ix.2")]
  AS_seg_cols = cbind( seg_colors[AS.seg.ix[,1]], seg_colors[AS.seg.ix[,2]] )

  if( plot.seg.sem ) 
  {
     AS_seg_sem = allele.segs[ , c("A1.Seg.sem", "A2.Seg.sem") ]
  }
  
  if( log2CR )
  {
     y.fn = function( y ) { log(y,2) }
  } else
  {
     y.fn = function( y ) { y }
  } 

  for (s in seq_len(nrow(allele.segs))) 
  {
    seg.crds <- as.numeric(c(allele.segs[s, "Start.bp"],
                             allele.segs[s, "End.bp"]))
    chr <- as.integer(allele.segs[s, "Chromosome"])
    genome.crds <- chr.offsets[chr] + seg.crds / chr.lens[chr] * chr.w[chr]
    means <- c(allele.segs[s, "A1.Seg.CN"],
               allele.segs[s, "A2.Seg.CN"], NA)
    
    het.seg.col.1 <- AS_seg_cols[s, 1]
    het.seg.col.2 <- AS_seg_cols[s, 2]

## make seg widths 95% CI for seg-mean
    min.w.d = 0.015/2  ## minimum possible width - make sure segs stay visible
    if( plot.seg.sem )
    {
       w.d_1 = max( AS_seg_sem[s,1] * 1.96, min.w.d )
       w.d_2 = max( AS_seg_sem[s,2] * 1.96, min.w.d )
    }
    else {  w.d_1 = w.d_2 = 0.025/2 }

    ybottom_1 =  means[1] - w.d_1
    ytop_1 = means[1] + w.d_1

    ybottom_2 = means[2] - w.d_2
    ytop_2 = means[2] + w.d_2

    if( !plot.total.CN )
    {
         rect(xleft = genome.crds[1], ybottom=y.fn(ybottom_1), xright = genome.crds[2],      ## minor CN
             ytop=y.fn(ytop_1) , border = NA, col = het.seg.col.1)

         rect(xleft = genome.crds[1], ybottom=y.fn(ybottom_2) , xright = genome.crds[2],     ## major CN
            ytop=y.fn(ytop_2), border = NA, col = het.seg.col.2)
    }
  }

  if( plot.total.CN )
  {
     for (s in seq_len(nrow(total.seg.dat))) 
     {
        seg.crds <- as.numeric(c(total.seg.dat[s, "Start.bp"], total.seg.dat[s, "End.bp"]))
        chr <- as.integer(total.seg.dat[s, "Chromosome"])
        genome.crds <- chr.offsets[chr] + seg.crds / chr.lens[chr] * chr.w[chr]
        tot.CR = total.seg.dat[s, "copy_num"]/2

        min.w.d = 0.015/2 
        if( plot.seg.sem )
        {
           w.d = max( total.seg.dat[s,"seg_sigma"] * 1.96, min.w.d )
        } else {  w.d = 0.025/2 }

        ybottom =  tot.CR - w.d
        ytop = tot.CR + w.d
 
        rect(xleft = genome.crds[1], ybottom=y.fn(ybottom), xright = genome.crds[2],      ## tot CN
             ytop=y.fn(ytop) , border = NA, col = tot.seg.colors[s] )
     }
  }

  if (plot.model.fit)
  {
     msg = get_SCNA_model_fit_plot_strings( mode.info )
     for( i in 1:length(msg)) { 
        mtext(msg[[i]], line = 2-i+1, cex = 0.75, adj = 0)
     }
  }

  if( plot.comb.fit )
  {
  ## ABS CN at each level
    for (q in seq_along(comb)) 
    {
       abline(h = y.fn(comb[q]), lwd = 0.5, lty = 3, col = comb.color)
    } 
  }
}


GetSegColsByAllelicBalance <- function( obs, allele.segs, seg.phase=NULL) 
{
  seg.col <- matrix(NA, nrow=nrow(allele.segs), ncol = 2)

  if (is.null(seg.phase)) 
  {
    f_vals = allele.segs[,"A1.Seg.CN"] / rowSums( allele.segs[,c("A1.Seg.CN", "A2.Seg.CN")] )
    seg.col =  499 * cbind(1-f_vals, f_vals) + 250
  }
  else
  {
     for (s in seq_len(nrow(allele.segs))) {
       e.state <- round((seg.phase[s, ] / 3) * 999, 0) + 1
       seg.col[s, 1] <- e.state[1]
       seg.col[s, 2] <- e.state[2]
     }
  }

  hom.color <- "darkgrey"
  het.1.color <- "red"
  het.2.color <- "blue"
  mid.color <- "darkviolet"
  pal <- colorRampPalette(c(hom.color, het.1.color, mid.color,
                            het.2.color, hom.color))
  cols <- pal(1000)

  seg_colors = c( cols[seg.col[, 1]], cols[seg.col[, 2]] )
#  AS.seg.cols = cbind( cols[seg.cols[, 1]], cols[seg.cols[, 2]] )
  
  return(seg_colors)
}


## Top-level called functions
PlotHscrAndSeghist <- function(allele.segs, seg_colors, max_CR, min_CR=-0.05, plot.genome=TRUE, plot.hist=FALSE, plot.abs.fit=FALSE, comb=NA, mode.info=NA, Wq0=NA, comb.ab=NA, fit.color=NA, plot.seg.sem=FALSE, plot.total.CN=FALSE, tot.seg.colors=NA, total.seg.dat=NA, y.lab=NA, log2CR=FALSE, label.chrs=TRUE, SID_label=NA) 
{
#  obs = seg.dat[["obs.scna"]]
#  allele.segs = get_hom_pairs_segtab( seg.dat )

  d.mar <- par("mar")
  ## eliminate right margins
  par(mar = c(d.mar[c(1:2)], d.mar[3], 0))

  if( plot.genome ) 
  {
     if( is.na(y.lab)) 
     { 
        y.lab = ifelse( plot.total.CN, "Copy ratio", "Allelic copy ratio") 
        if( log2CR ) { y.lab = paste( "Log2 ", y.lab, sep="") }
     }

     GenomeHscrSegPlot(allele.segs, seg_colors, y.lab=y.lab, y.min=min_CR, y.max=max_CR, plot.model.fit=!is.na(mode.info), plot.comb.fit=plot.abs.fit, mode.info=mode.info, comb=comb, comb.color=fit.color, plot.seg.sem=plot.seg.sem, plot.total.CN=plot.total.CN, total.seg.dat=total.seg.dat, tot.seg.colors=tot.seg.colors, log2CR=log2CR, label.chrs=label.chrs ) 
   
     if( !is.na(SID_label)) { 
        mtext( SID_label, side=3, line=0.5, adj=0, cex=par("cex")+0.1, font=2 )
     }
  }
  
  ##  seg hist plots
  bin.w = 0.04  ## for histograms
  ## no left margin
  par(mar = c(d.mar[1], 0, d.mar[3], d.mar[4]))
  
  if (plot.hist) 
  {
    if( !plot.total.CN )
    {
       seg_means = c(allele.segs[,"A1.Seg.CN"], allele.segs[,"A2.Seg.CN"])
       seg_W = c( allele.segs[,"W"], allele.segs[,"W"] )
    }
    else 
    {
       seg_means = total.seg.dat[,"copy_num"]/2
       seg_W = total.seg.dat[,"W"]
       comb.ab=NA
       seg_colors = tot.seg.colors
    }
 
    if( log2CR ) { seg_means = log(seg_means, 2) }
    if( log2CR ) { comb = log(comb, 2) }

    plot_ABS_seg_hist(seg_means, seg_W, copy_ratio_label="", seg_colors, min_CR, max_CR, bin.w, sideways=TRUE )

#    mtext(text = "Fraction of reference genome", line=par("mgp")[1], side = 1, cex=par("cex.axis")*par("cex") )
#    abline(h = 0, lty = 3, lwd = 2)
#    mtext("Summary histogram", side = 3, line = 0, adj = 0, cex = par("cex.axis"))

    if (plot.abs.fit) 
    {
       plot_ABS_comb_fit( sideways=TRUE, col=fit.color, max_CR=max_CR, plot.model.fit=FALSE, comb=comb, mode.info=mode.info, Wq0=Wq0, comb.ab=comb.ab )
    }

    if( !plot.genome & !is.na(SID_label)) { 
        mtext( SID_label, side=3, line=0.5, adj=0, cex=par("cex")+0.1, font=2 )
    }

  } 

## reset plot margins
   par(mar=d.mar)
}



