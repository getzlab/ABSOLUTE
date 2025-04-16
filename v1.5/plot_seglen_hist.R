PlotSeglenHist <- function(D, W, seg_colors=NA, color.by=NA, color.range=NA, x.max=NA, x.min=NA,
                             y.max=NA, bin.w=0.05, order.by=NA, use.pal=NA,
                             data.whiskers=TRUE,
                             xlab="Observed allelic copy-ratio", x.axis=TRUE,
                             y.axis=TRUE, x.ax.labs=TRUE, y.ax.labs=TRUE,
                             xlim=NA, ylab="Genomic fraction", border=TRUE,
                             sideways=FALSE,add=FALSE) {
  
  # Sometimes when the tumor has no copy-number alterations, ABSOLUTE attempts to plot a chart with no data.
  # To prevent a crash, I just use an arbitrary color.
  if( all(is.na(seg_colors)) & is.na(color.by)) { seg_colors = c(1) }

 
  P <- length(D)
  
  if (is.na(x.max)) {
    x.max <- max(D)
  }
 
  if( is.na(x.min)) { 
    x.min = min(D)
  }
  
  if (length(bin.w) == 1) {
    res <- hist(D, breaks=seq(x.min, x.max + bin.w, bin.w), plot=FALSE)
    breaks <- res[["breaks"]]
  }  else {
    breaks <- bin.w
  }
  
  heights <- matrix(0, ncol=length(breaks), nrow=length(D))
  seg.colors <- matrix(0, ncol=length(breaks), nrow=length(D))
  
  if (is.na(use.pal)) {
    col.scale <- 1000
    use.pal <- heat.colors(col.scale)
  } else {
    col.scale <- length(use.pal)
  }
  
  if( !is.na(color.by) )
  {
     if (is.na(color.range)) {
       color.range <- range(color.by, na.rm=TRUE)
     }
     pal.idx <- floor((color.by - color.range[1]) /
                      (color.range[2]-color.range[1]) * (col.scale-1))
     
    for (p in seq_len(P)) {
      bin <- max(which(breaks <= D[p]))
      heights[p, bin] <- W[p]
      seg.colors[p, bin] <- use.pal[pal.idx[p] + 1]
    }
  }

  if( !is.na(seg_colors) ) ## input arg
  {

     for (p in seq_len(P)) {
         bin <- max(which(breaks <= D[p]))
         heights[p, bin] <- W[p]
         seg.colors[p, bin] <- seg_colors[p]
       }
  }
  
  ## un-modeled segments
  seg.colors[is.na(seg.colors)] <- "grey" 
  
  ## sort cols 
  for (i in seq_len(ncol(heights))) {
    if (is.na(order.by)) {
      res <- sort(heights[, i], decreasing=TRUE, index.return=TRUE)
    } else { 
      res <- sort(order.by, index.return=TRUE)
    }   
    
    heights[, i] <- heights[res[["ix"]], i]
    seg.colors[, i] <- seg.colors[res[["ix"]], i]
  }
  
  seg.colors[heights == 0] <- NA
  heights[heights == 0] <- NA
  
  MyBarplot(heights, seg.colors, ylab, xlab, x.max, y.max,
             add=add, breaks=breaks, xlim=xlim, border=border,
             sideways=sideways)
  
  ## show data whiskers
  if (data.whiskers) {
    side <- ifelse(!sideways, 1, 2)
    axis(side=side, at=D, labels=FALSE, line= -0.5, col="grey")
    axis(side=side, at=D[!is.na(color.by)], labels=FALSE, line= -0.5, col="red")
    x.axis.line <- 1
  } else{
    x.axis.line <- NA
  }
  
  ## numbered data axis
  if (x.axis) {
    if (!sideways) {
      axis(side=1, line=x.axis.line, labels=x.ax.labs)
    } else {
      axis(side=2, line=x.axis.line, labels=x.ax.labs)
    }
  }
  
  if (y.axis) {
    if (!sideways) {
      par("ylog"=FALSE)
      axis(side=2, labels=y.ax.labs, las=2)
      par("ylog"=FALSE)
    } else {
      par("xlog"=FALSE)
      axis(side=1, labels=y.ax.labs, las=1)
      par("xlog"=FALSE)
    }
  }  
}

MyBarplot <- function(height.mat, col.mat, ylab, xlab, max.data, max.y,
                      add=FALSE, break.sz=0.05, breaks, xlim=NA, ylim=NA,
                      border=TRUE, sideways=FALSE) {
  cur.stack.height <- rep(0, ncol(height.mat))
  stack.heights <- apply(height.mat, 2, sum, na.rm=TRUE)
  max.height <- max(stack.heights, na.rm=TRUE)
  if (!is.na(max.y)) {
    max.height <- max.y
  }
  
  if (!add) {
    if (is.na(xlim)) {
      xlim <- c(0,max.data+break.sz / 2)
    }

    if (!sideways) {
      plot(x=0,y=0, xlim=xlim, ylim=c(0,max.height), type="n", xlab=xlab,
           ylab=ylab, bty="n", axes=FALSE)
      abline(h=0)
    } else {
      plot(x=0,y=0, ylim=xlim, xlim=c(0,max.height), type="n", xlab=xlab,
           ylab=ylab, bty="n", axes=FALSE)
      abline(v=0)
    }
  }

  for (i in seq_len(nrow(height.mat))) {
    for (j in seq_len(ncol(height.mat))) {
      if (!is.na(height.mat[i, j])) {
        ## omit border for tiny segments
        if ((border == TRUE) & height.mat[i, j] / max.height > 0.01) {
          border.color <- 1
        } else { 
          border.color <- col.mat[i, j] 
        }
        if (!sideways) {
          rect(xleft=breaks[j], ybottom=cur.stack.height[j],
               xright=breaks[j+1], ytop=cur.stack.height[j]+height.mat[i,j],
               col=col.mat[i,j], border=border.color, lwd=0.5)
        } else {
          rect(xleft= cur.stack.height[j], 
               ybottom=breaks[j], 
               xright=cur.stack.height[j]+height.mat[i, j],
               ytop=breaks[j+1], 
               col=col.mat[i, j], border=border.color, lwd=0.5)
        }
        cur.stack.height[j] <- cur.stack.height[j]+height.mat[i, j]
      } 
    }
  }
}

