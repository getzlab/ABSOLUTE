## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GetMutBetaDensities <- function(mut.dat, n.grid=100) 
{
  mut.grid = matrix(0, nrow=nrow(mut.dat), ncol=n.grid)
  
  cov = mut.dat[, "alt"] + mut.dat[, "ref"]
  af = mut.dat[, "alt"] / cov
  
  grid.vals = seq_len(n.grid) / (n.grid + 1)
  
  for (i in seq_len(nrow(mut.grid))) {
    mut.grid[i, ] = dbeta(grid.vals, cov[i] * af[i] + 1, 
                          cov[i] * (1 - af[i]) + 1) * 1 / n.grid
  }

  bad_rows = which( apply( is.nan(mut.grid), 1, sum ) > 0 )
  if(length(bad_rows)>0 ) { stop() }
  
  return(mut.grid)
}

DrawMutBetaDensities <- function(beta.grid, pr.clonal, hz.del.flag, cols,
                                 draw.indv=TRUE, draw.total=TRUE) {
  n.grid <- ncol(beta.grid)
  grid.vals <- seq_len(n.grid) / (n.grid + 1)
  pr.subclonal = 1-pr.clonal
  pr.subclonal[pr.subclonal<0] = 0 ## round-off error
  
  pal <- colorRampPalette(cols)
  colpal <- pal(1000)
  col.scale <- length(colpal)
  color.range <- c(0, 1)
  color.by <- pr.clonal
  pal.idx <- floor((color.by - color.range[1]) / (color.range[2] - color.range[1]) * 
    (col.scale - 1)) + 1
  mut.colors <- colpal[pal.idx]
  
  mut.colors[hz.del.flag] <- "navy"
  
  if (draw.indv) {
    for (i in seq_len(nrow(beta.grid))) {
      lines(grid.vals, beta.grid[i, ], col=mut.colors[i])
    }
  }
  
  ## Pr-weighted 
  clonal.grid <- matrix(0, nrow=nrow(beta.grid), ncol=ncol(beta.grid))
  sc.grid <- clonal.grid
  for (i in seq_len(nrow(beta.grid))) {
    clonal.grid[i, ] <- beta.grid[i, ] * pr.clonal[i]
    sc.grid[i, ] <- beta.grid[i, ] * pr.subclonal[i]
  }
  
  if (draw.total) {
    lines(grid.vals, colSums(clonal.grid)/max(colSums(clonal.grid)), col=cols[2], lty=4)
    lines(grid.vals, colSums(sc.grid)/max(colSums(sc.grid)), col=cols[1], lty=4)
  }
}

Mut_AF_plot = function(mut.dat, SSNV_cols, mode.color, draw.indv)
{
  cov <- mut.dat[, "alt"] + mut.dat[, "ref"]

  pr.clonal <- mut.dat[, "Pr_somatic_clonal"]
  pr.cryptic.SCNA <- mut.dat[, "Pr_cryptic_SCNA"]
  SSNV_skew <- mut.dat[1, "SSNV_skew"]

  n_grid = 300
  af_post_pr = GetMutBetaDensities(mut.dat, n_grid)
  grid = 1:n_grid / (n_grid + 1)
  grid_mat = matrix(grid, nrow=nrow(af_post_pr), ncol=ncol(af_post_pr), byrow=TRUE)
  
  plot(0, type="n", bty="n", main="", xlab="Fraction of alternate reads", ylab="Density", 
       xlim=c(0, 1), ylim=c(0, 1), las=1 ) #max(colSums(af_post_pr))
  hz.del.flag <- mut.dat[, "q_hat"] == 0
  DrawMutBetaDensities(af_post_pr, pr.clonal, hz.del.flag, SSNV_cols, draw.total=TRUE,
                       draw.indv=draw.indv)

#  if(!is.na(SSNV_skew)){
#     msg = paste("skew = ", round(SSNV_skew,3), sep="")
#     mtext( text=msg, side=3, adj=1  )
#  }

  alpha <- mut.dat[1, "purity"]
  abline(v=alpha / 2, lwd=0.5, lty=3, col=mode.color)
  msg <- expression(hat(alpha) / 2)
  mtext(text=msg, side=3, at=alpha/2, col=mode.color, line=0.2, cex=par("cex") * par("cex.axis"))

 
  abline(v=alpha*SSNV_skew/2, lwd=0.5, lty=3, col="black")
  msg <- expression(hat(f_s) * hat(alpha) / 2)
  mtext(text=msg, side=3, at=alpha*SSNV_skew/2, col="black", line=0.2, cex=par("cex") * par("cex.axis"))

  return( list("af_post_pr"=af_post_pr, "grid_mat"=grid_mat))
}      



draw_mut_multiplicity_densities = function(mut_pr, grid, pr_clonal, pr_cryptic_SCNA, x_lim, 
                                   xlab, cols, draw_indv=TRUE, draw_total=TRUE, 
                                   add=FALSE, y_lim=NA) 
{  
  get_grid_combined_mut_densities = function(mut_pr, pr_clonal, grid, x_lim) 
  {
    bin_w = x_lim / 100
    breaks=seq(0, x_lim, by=bin_w)
    mult_grid = breaks
  
  ## pr-weighted 
    grid_dens = matrix(0, nrow=nrow(mut_pr), ncol=length(mult_grid))
  
    for (i in seq_len(nrow(mut_pr))) 
    {
      x = grid[i, ] 
#      y = mut_pr[i, ] * pr_clonal[i] 
      y = mut_pr[i, ]
      if (sum(!is.na(y)) > 2) {
        grid_dens[i, ] = approx(x, y, xout=mult_grid)$y
      }
    }

    y_lim = colSums(grid_dens, na.rm=TRUE)
    return(list(grid_dens=grid_dens, MULT_GRID=mult_grid, YLIM=y_lim))  
  }

  pr_subclonal = 1 - pr_clonal
  pr_subclonal[pr_subclonal < 0] = 0  ## round-off error

  pal = colorRampPalette(cols)
  colpal = pal(1000)
  
  col_scale = length(colpal)
  color_range = c(0, 1)
  color_by = pr_clonal
  pal_idx = floor((color_by - color_range[1]) / (color_range[2] - color_range[1]) * (col_scale-1)) + 1
  mut_colors = colpal[pal_idx]
  
  ix = pr_cryptic_SCNA > 0.5
  mut_colors[ix] = "mediumorchid2"
  
  # stack up variable-grid rescaled densities
  res = get_grid_combined_mut_densities(mut_pr, pr_clonal, grid, x_lim)
  clonal_dens = res$grid_dens

  res = get_grid_combined_mut_densities(mut_pr, pr_subclonal, grid, x_lim)
  sc_dens = res$grid_dens 
  mult_grid= res$MULT_GRID
  
  if (!add) {
    plot( 0, type="n", bty="n", main="", xlab=xlab, ylab="Density", xlim=c(0,x_lim), ylim=c(0, y_lim), las=1)
  }
  
  if (draw_indv) {
    for (i in seq_len(nrow(mut_pr))) {
#      lines(grid[i, ], mut_pr[i, ], col=mut_colors[i])
      lines(mult_grid, clonal_dens[i, ], col=mut_colors[i])
    }
  }
  
  clonal_dens[is.na(clonal_dens)] = 0
  sc_dens[is.na(sc_dens)] = 0
  
  if (draw_total) {
  
	##edit
	#print(colSums(clonal_dens*pr_clonal))
	ncl=colSums(clonal_dens*pr_clonal)
	nsbcl=colSums(sc_dens*pr_subclonal)

    lines(mult_grid, ncl/max(ncl), col=cols[2], lty=4)
    lines(mult_grid, nsbcl/max(nsbcl), col=cols[1], lty=4)

#	save(ncl, nsbcl, file = "clon_test.RData")


	}
	

}


multiplicity_plot = function( seg.dat, mut.dat, af_post_pr, grid_mat, SSNV_cols, mode.color, draw.indv, verbose=FALSE )
{
  hz.del.ix <- mut.dat[, "q_hat"] == 0
  SC_CN.ix = !is.na(mut.dat[,"H1"]) ## SSNVs on subclonal SCNAs
  nix = hz.del.ix | SC_CN.ix |  mut.dat[,"alt"]==0  # drop force-called muts

  mut.dat <- mut.dat[!nix, , drop=FALSE]
  af_post_pr = af_post_pr[!nix,, drop=FALSE]
  grid_mat = grid_mat[!nix,, drop=FALSE]

  if (nrow(mut.dat) == 0) {
    if (verbose) { print("No valid SSNVs to plot") }
    frame()
    return()
  }

 
  res = get_SSNV_on_clonal_CN_multiplicity_densities( seg.dat, mut.dat, af_post_pr, grid_mat, verbose=verbose )


  mult_dens = res[["mult_dens"]]
  mult_grid = res[["mult_grid"]]
 
  pr.clonal <- mut.dat[, "Pr_somatic_clonal"]
  pr.cryptic.SCNA <- mut.dat[, "Pr_cryptic_SCNA"]
  SSNV_skew <- mut.dat[1, "SSNV_skew"]


  mult.xlim <- 2.5
  draw_mut_multiplicity_densities(mult_dens, mult_grid, pr.clonal, pr.cryptic.SCNA, x_lim=mult.xlim,
                          cols=SSNV_cols, xlab="SSNV multiplicity", draw_total=TRUE,
                          draw_indv=draw.indv, y_lim=1 )  ##log(max(colSums(mult_dens))) changed to max=1

  abline( v=1, lty=3, lwd=0.5, col=mode.color )  ## highlight integer multiplicities
  abline( v=2, lty=3, lwd=0.5, col=mode.color )


  abline( v=SSNV_skew*1, lty=3, lwd=0.5, col=1 )  ## highlight skew at integer multiplicities
  abline( v=SSNV_skew*2, lty=3, lwd=0.5, col=1 )

}

draw_grid_mut_densities = function(mut_pr, grid, col, mut_colors, scale_total, draw_total=TRUE, draw_indv=TRUE, add=FALSE, y_lim=NA, lty=1) 
{  
  if (!add) {
    plot( 0, type="n", bty="n", main="", xlab="CCF", ylab="Density", xlim=c(0,1), ylim=c(0, y_lim), las=1)
  }
  
  if (draw_indv) {
    for (i in seq_len(nrow(mut_pr))) {
      lines(grid, mut_pr[i, ], col=mut_colors[i], lty=lty)
    }
  }
  
  if (draw_total) {
    combined_dens = colMeans(mut_pr) + scale_total
    lines(grid, combined_dens, col=col, lty=4)
  }
}

get_mut_cols = function( cols, color_by, SSNV.amp.ix )
{
  pal = colorRampPalette(cols[c(1,2)] )
  colpal = pal(1000)
  col_scale = length(colpal)
  color_range = c(0, 1)
  pal_idx = floor((color_by - color_range[1]) / (color_range[2] - color_range[1]) * (col_scale-1)) + 1
  mut_colors = colpal[pal_idx]

  mut_colors[ SSNV.amp.ix ] = cols[3]

  return(mut_colors)
}

get_SSNV_plot_colors = function( mut.dat, SSNV_cols )
{
  H123.SSNV.ix = !is.na(mut.dat[,"H1"]) ## SSNVs on subclonal SCNAs
  SSNV_pr_clonal <- mut.dat[, "Pr_somatic_clonal"]

  SSNV.amp.ix = mut.dat[,"modal_q_s"] > 1

  SSNV_mut_cols = get_mut_cols( SSNV_cols, SSNV_pr_clonal[!H123.SSNV.ix], SSNV.amp.ix )

#  H123.SSNV_mut_cols = get_mut_cols( H123_SSNV_cols, SSNV_pr_clonal[H123.SSNV.ix] )
  H123.SSNV_mut_cols = rgb( "red"=mut.dat[H123.SSNV.ix,"H1"], "green"=mut.dat[H123.SSNV.ix,"H2"], "blue"=mut.dat[H123.SSNV.ix,"H3"] )


  cols = rep(NA, nrow(mut.dat) )
  cols[ H123.SSNV.ix ] = H123.SSNV_mut_cols
  cols[ !H123.SSNV.ix ] = SSNV_mut_cols

  return(cols)
}


SSNV_CCF_plot = function( seg.dat, mut.dat, SSNV_ccf_dens, mode.ix, SSNV_cols, draw.indv )
{
  H123.SSNV.ix = !is.na(mut.dat[,"H1"]) ## SSNVs on subclonal SCNAs
  SSNV_pr_clonal <- mut.dat[, "Pr_somatic_clonal"]
  ccf_grid = as.numeric(colnames(SSNV_ccf_dens))
  H123_SSNV_cols = c("darkslateblue", "coral", "seagreen3")

  clonal_SCNA_SSNV_ccf_dens = SSNV_ccf_dens[ !H123.SSNV.ix, ,drop=FALSE ]
  H123.SSNV_ccf_dens = (SSNV_ccf_dens)[ H123.SSNV.ix, , drop=FALSE ]

  clonal_SCNA_SSNV_ymax = max(clonal_SCNA_SSNV_ccf_dens, na.rm=TRUE)
  H123.SSNV_ymax = max( (H123.SSNV_ccf_dens))

  y_lim = max( clonal_SCNA_SSNV_ymax, H123.SSNV_ymax, na.rm=TRUE  )
  scale_total = y_lim/2

  SSNV.amp.ix = mut.dat[,"modal_q_s"] > 1
  SSNV_mut_cols = get_mut_cols( SSNV_cols, SSNV_pr_clonal[!H123.SSNV.ix], SSNV.amp.ix )
  H123.SSNV_mut_cols = get_mut_cols( H123_SSNV_cols, SSNV_pr_clonal[H123.SSNV.ix], SSNV.amp.ix )

## SSNVs on clonal and subclonal SCNAs
  draw_grid_mut_densities( clonal_SCNA_SSNV_ccf_dens, ccf_grid, col=SSNV_cols[2], mut_colors=SSNV_mut_cols, scale_total, draw_total=TRUE, draw_indv=TRUE, add=FALSE, y_lim=y_lim )

  draw_grid_mut_densities( H123.SSNV_ccf_dens, ccf_grid, col=H123_SSNV_cols[2], mut_colors=H123.SSNV_mut_cols, scale_total, draw_total=TRUE, draw_indv=TRUE, add=TRUE, y_lim=y_lim )


  legend( x='topleft', legend=c("Clonal SSNVs", "Subclonal SSNVs", "Clonal SSNVs on subclonal SCNAs", "Subclonal SSNVs on subclonal SCNAs"), col=c( rev(SSNV_cols), rev(H123_SSNV_cols) ), lty=1, lwd=1.5, bty="n" )
}


SSNV_on_subclonal_SCNA_CCF_plot = function( seg.dat, mut.dat, SSNV_ccf_dens, mode.ix, SSNV_cols, draw.indv )
{
  H123.SSNV.ix = !is.na(mut.dat[,"H1"]) ## SSNVs on subclonal SCNAs
  if( sum(H123.SSNV.ix) == 0 ) { frame(); return() }

  mut.dat = mut.dat[H123.SSNV.ix,,drop=FALSE]
  SSNV_ccf_dens = SSNV_ccf_dens[H123.SSNV.ix, ,drop=FALSE ]

  SC_Aq_a = mut.dat[,"SC_Aq_a"]
  SC_Aq_d = mut.dat[,"SC_Aq_d"]
  loss.ix = SC_Aq_d < SC_Aq_a

  ccf_grid = as.numeric(colnames(SSNV_ccf_dens))
  H_cols= c(2,3,4)
  y_lim = max( SSNV_ccf_dens, na.rm=TRUE  )
  scale_total = y_lim/2
  mut_colors = rgb( "red"=mut.dat[,"H1"], "green"=mut.dat[,"H2"], "blue"=mut.dat[,"H3"] )

  draw_grid_mut_densities( SSNV_ccf_dens, ccf_grid, col=NA, mut_colors=mut_colors, scale_total, draw_total=FALSE, draw_indv=TRUE, add=FALSE, y_lim=y_lim )

  draw_grid_mut_densities( SSNV_ccf_dens[loss.ix,,drop=FALSE], ccf_grid, col=NA, mut_colors=rep("black", sum(loss.ix)), scale_total, draw_total=FALSE, draw_indv=TRUE, add=TRUE, y_lim=y_lim, lty=3 )

  mtext( text="SSNVs on subclonal SCNAs", side=3, line=0, adj=0, cex=0.7)
  legend( x='topleft', legend=c("H1: Ancestral SSNV in trans", "H2: Ancestral SSNV in cis", "H3: Derived SSNV", "Subclonal deletion SCNA"), col=c(H_cols,"black"), lty=c(1,1,1,3), lwd=1.5, bty="n" )
}

get_grid_dens_95CI = function( dens, grid )
{
  mult_dens = dens
  mult_grid = grid

# calc mult mode and 95% CI
   mult_dens= mult_dens/ rowSums(mult_dens)
   mult_hat = rep(NA, nrow(dens))
   mult_ci95 = matrix( NA, nrow=nrow(dens), ncol=2)

   for( i in 1:nrow(dens))
   {
      if( all(is.na(mult_dens[i,]))) { next }

      mult_hat[i] = mult_grid[i, which.max(mult_dens[i,]) ]

#      if( is.na(mult_hat[i])) { next }

      ecdf = cumsum(mult_dens[i,] )
      if(ecdf[1]==1) { mult_ci95[i,]=c(0,0.1) }
      else
      {
         mult_ci95[i, ] = approx(ecdf, y=mult_grid[i,], xout=c(0.025, 0.975))$y
      }
   }

   return( cbind(mult_ci95, mult_hat) )
}


## External called function
plot_SSNVs_on_genome = function( SSNV_model, SSNV_colors, mut.dat, seg.dat, i, mode.color, min.cov=3, verbose=FALSE )
{
   plot_SSNV_genome_vs_CI95 = function( mut.dat, colors, CI95 )
   {
      chrpos = cbind( mut.dat[,"Chromosome"], mut.dat[,"Start_position"] )
      genome_coords = get_genome_coords( chrpos )
      ycrd = CI95[,3]  ## modal point est
      points( genome_coords, ycrd, pch=16, col=colors, cex=0.65 )

      segments( x0=genome_coords, y0=CI95[,1], y1=CI95[,2], col=colors, lwd=0.5 )
   }

   nix1 = mut.dat[,"alt"]==0  # drop force-called muts

   SC_CN.ix = !is.na(mut.dat[,"H1"]) ## SSNVs on subclonal SCNAs
   nix1.SC = mut.dat[ SC_CN.ix, "alt"]==0
   H.123.ssnv.ccf.dens = SSNV_model[["SSNV_on_subclonal_SCNA_res"]][["H.123.ssnv.ccf.dens"]] [, !nix1.SC,,drop=FALSE]
   mut.dat = mut.dat[!nix1,, drop=FALSE]
   SC_CN.ix = !is.na(mut.dat[,"H1"]) ## SSNVs on subclonal SCNAs



   n_grid=100
   af_post_pr = GetMutBetaDensities(mut.dat, n_grid)
   grid = 1:n_grid / (n_grid + 1)
   grid_mat = matrix(grid, nrow=nrow(af_post_pr), ncol=ncol(af_post_pr), byrow=TRUE)

   if( any(!SC_CN.ix)) 
   {
      res = get_SSNV_on_clonal_CN_multiplicity_densities( seg.dat, mut.dat[!SC_CN.ix,,drop=FALSE], af_post_pr[!SC_CN.ix,,drop=FALSE], grid_mat[!SC_CN.ix,,drop=FALSE], verbose=verbose )
      C_mult_dens = res[["mult_dens"]]
      C_mult_grid = res[["mult_grid"]] 
      C_mult_ci95 = get_grid_dens_95CI( C_mult_dens, C_mult_grid )

      colors = get_SSNV_plot_colors( mut.dat[!SC_CN.ix,,drop=FALSE], SSNV_colors )
      plot_SSNV_genome_vs_CI95( mut.dat[!SC_CN.ix,,drop=FALSE], colors, C_mult_ci95 )
   }

   SC_SSNV_colors = c("red", "green", "blue", "black")

   if( any( SC_CN.ix ))
   {
      mut.H.Prs = mut.dat[ SC_CN.ix, c("H1", "H2", "H3", "H4"), drop=FALSE ]
      grid_mat = matrix( SSNV_model[["ccf_grid"]], nrow=sum(SC_CN.ix), ncol=length(SSNV_model[["ccf_grid"]]), byrow=TRUE)

      for( i in 1:(dim(H.123.ssnv.ccf.dens)[1]) )
      {
         CCF_dens = matrix( H.123.ssnv.ccf.dens[i,,], nrow=sum(SC_CN.ix) )
         SC_mult_ci95 = get_grid_dens_95CI( CCF_dens, grid_mat )
         
         r = col2rgb(SC_SSNV_colors[i])[1]
         g = col2rgb(SC_SSNV_colors[i])[2]
         b = col2rgb(SC_SSNV_colors[i])[3]

         TRANS = 255 * mut.H.Prs[,i]
         colors = rgb( r, g, b, TRANS, maxColorValue=255 )

         plot_SSNV_genome_vs_CI95( mut.dat[SC_CN.ix,,drop=FALSE], colors, SC_mult_ci95 )

      }
   }
}




## External called function
PlotSomaticMutDensities <- function(mut.dat, seg.dat, mode.ix, mode.color, min.cov=1, verbose=FALSE)
{
  cov <- mut.dat[, "alt"] + mut.dat[, "ref"]
  ix <- cov > min.cov
  mut.dat <- mut.dat[ix, , drop=FALSE]

  SSNV_ccf_dens = seg.dat[["mode.res"]][["SSNV.ccf.dens"]][mode.ix, , ]
  SSNV_ccf_dens = SSNV_ccf_dens[ix,,drop=FALSE]

  if (nrow(mut.dat) == 0) {
    if (verbose) {
      print("No valid SSNVs to plot")
    }
    for( i in 1:4){ frame() }
    return()
  }

  SSNV_cols = c("dodgerblue", "darkgrey", "seagreen3")   ## SC, clonal, mult>1

  draw.indv <- nrow(mut.dat) < 500

## 1st plot: mut VAFs with clonal / SC colors and purity line
  res = Mut_AF_plot(mut.dat, SSNV_cols, mode.color, draw.indv)

## 2nd plot: Plot rescaling to multiplicity
  multiplicity_plot( seg.dat, mut.dat, res[["af_post_pr"]], res[["grid_mat"]], SSNV_cols, mode.color,  draw.indv, verbose )

## 3rd plot: combined SCNAs and muts
  SSNV_CCF_plot( seg.dat, mut.dat, SSNV_ccf_dens, mode.ix, SSNV_cols, draw.indv )

## 4th plot: CCF of muts on subclonal SCNAs by H* class
  SSNV_on_subclonal_SCNA_CCF_plot( seg.dat, mut.dat, SSNV_ccf_dens, mode.ix, SSNV_cols, draw.indv )
}



