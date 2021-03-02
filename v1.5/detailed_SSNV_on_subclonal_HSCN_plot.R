detailed_SSNV_on_subclonal_HSCN_plot = function( SSNV_model, seg.dat, mut.dat, mode.ix, min.cov=3, verbose=FALSE )
{
   SSNV_on_subclonal_SCNA_res = SSNV_model[["SSNV_on_subclonal_SCNA_res"]]

#  SSNV_cols = c("maroon", "olivedrab4") 
   SSNV_ccf_dens = seg.dat[["mode.res"]][["SSNV.ccf.dens"]][mode.ix, , ]

#  cov <- mut.dat[, "alt"] + mut.dat[, "ref"]
#  ix <- cov > min.cov
#  cov = cov[ix]
#  mut.dat <- mut.dat[ix, , drop=FALSE]
#  SSNV_ccf_dens = SSNV_ccf_dens[ix,,drop=FALSE]

  if (nrow(mut.dat) == 0) {
    if (verbose) {
      print("No valid SSNVs to plot")
    }
    return()
  }


  sc_tab = seg.dat$mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[mode.ix, , ]
  SCNA_pr_subclonal = sc_tab[, "Pr_subclonal"]
  SCNA_ccf_dens = exp( seg.dat$mode.res$subclonal_SCNA_res$log_CCF_dens[mode.ix, , ] )



  H123.SSNV.ix = !is.na(mut.dat[,"H1"]) ## SSNVs on subclonal SCNAs
  if( sum(H123.SSNV.ix) == 0 ) { frame(); return() }

  mut.dat = mut.dat[H123.SSNV.ix,,drop=FALSE]
  SSNV_ccf_dens = SSNV_ccf_dens[H123.SSNV.ix, ,drop=FALSE ]

  SC_Aq_a = mut.dat[,"SC_Aq_a"]
  SC_Aq_d = mut.dat[,"SC_Aq_d"]
  loss.ix = SC_Aq_d < SC_Aq_a

  ccf_grid = as.numeric(colnames(SSNV_ccf_dens))
  H_cols= c(2,3,4)
  mut_colors = rgb( "red"=mut.dat[,"H1"], "green"=mut.dat[,"H2"], "blue"=mut.dat[,"H3"] )



  H.123.ssnv.ccf.dens = SSNV_on_subclonal_SCNA_res[["H.123.ssnv.ccf.dens"]]
  CCF_dens = array( NA, dim=c(4, nrow(mut.dat), length(ccf_grid), length(ccf_grid) ) )
  for( i in 1:nrow(mut.dat) )
  {
    SC.ix = mut.dat[i,"SC.ix"]
    CCF_dens[1,i,,] = SCNA_ccf_dens[SC.ix,] %o% SSNV_ccf_dens[i,]

    CCF_dens[2,i,,] = SCNA_ccf_dens[SC.ix,] %o% H.123.ssnv.ccf.dens[1,i,]
    CCF_dens[3,i,,] = SCNA_ccf_dens[SC.ix,] %o% H.123.ssnv.ccf.dens[2,i,]
    CCF_dens[4,i,,] = SCNA_ccf_dens[SC.ix,] %o% H.123.ssnv.ccf.dens[3,i,]
  }

  par(mfrow=c(6,6))
  SCNA_SSNV_joint_posterior_CCF_plots( mut.dat, ccf_grid, CCF_dens, mut_colors )


#  draw_grid_mut_densities( SSNV_ccf_dens, ccf_grid, col=NA, mut_colors=mut_colors, scale_total, draw_total=FALSE, draw_indv=TRUE, add=FALSE, y_lim=y_lim )

#  draw_grid_mut_densities( SSNV_ccf_dens[loss.ix,,drop=FALSE], ccf_grid, col=NA, mut_colors=rep("black", sum(loss.ix)), scale_total, draw_total=FALSE, draw_indv=TRUE, add=TRUE, y_lim=y_lim, lty=3 )

#  mtext( text="SSNVs on subclonal SCNAs", side=3, line=0, adj=0, cex=0.7)
#  legend( x='topleft', legend=c("H1: Ancestral SSNV in trans", "H2: Ancestral SSNV in cis", "H3: Derived SSNV", "Subclonal deletion SCNA"), col=c(H_cols,"black"), lty=c(1,1,1,3), lwd=1.5, bty="n" )
}





## plot joint SCNA, SSVN CCF dens for each H..., and combined.

SCNA_SSNV_joint_posterior_CCF_plots = function( mut.dat, GRID, CCF_dens, mut_cols )
{
  source("~/CGA/R/Phylogic/pair_CCF_compare.R")


   mut_contour_plot = function( dens, GRID, level_scales, col, lwd=0.5 )
   {
      mode_dens = max(dens)
      CL = contourLines( x=GRID, y=GRID, z=dens, levels=level_scales * mode_dens)
   
      for( j in 1:length(CL) )
      {
         lines( CL[[j]]$x, CL[[j]]$y, col=col, lwd=0.5 )
      }
   }

   dmar = par("mar")
   nmar=dmar
   nmar[c(1,3)] = c(0,1.35)  # no bot, small top margin
   nmar[c(2,4)] = 0  # no left or right margin
   par(mar=nmar)
   
   level_scales = c( 0.99, 0.95, 0.8, 0.5, 0.1 )

   cov= mut.dat[,"alt"] + mut.dat[,"ref"]
   
   for( i in 1:nrow(mut.dat) )
   {
      plot( 0, type="n", bty="n", axes=FALSE, main="", ylab="", xlab="", xlim=c(0,1), ylim=c(0,1) )
      polygon( c(0,1,1,0), c(0,0,1,1), border=FALSE, col="grey95", add=TRUE )

#      mut_contour_plot( CCF_dens[2,i,,], GRID, level_scales, 2 )
#      mut_contour_plot( CCF_dens[3,i,,], GRID, level_scales, 3 )
#      mut_contour_plot( CCF_dens[4,i,,], GRID, level_scales, 4 )

      plot_trans_contours_2d_interval( CCF_dens[2,i,,], GRID, 2, level_scales )
      plot_trans_contours_2d_interval( CCF_dens[3,i,,], GRID, 3, level_scales )
      plot_trans_contours_2d_interval( CCF_dens[4,i,,], GRID, 4, level_scales )
      mut_contour_plot( CCF_dens[1,i,,], GRID, level_scales, mut_cols[i], lwd=1 )


 # mut.dat[i,"Variant_Classification"]

      mut_lab_1 = paste( mut.dat[i,"Hugo_Symbol"], " (", mut.dat[i,"alt"], "/", cov[i], "), ",  sep="" )
      mut_lab_2 = paste( "Aqd=", mut.dat[i,"SC_Aq_d"], ", Aqa=", mut.dat[i,"SC_Aq_a"], ", CAq=", mut.dat[i,"C_Aq"], sep="")

#      var_info = ifelse( mut.dat[i,"Variant_Classification"] %in% c("Missense_Mutation", "Nonsense_Mutation"), as.character(mut.dat[i,"Protein_Change"]),  as.character(mut.dat[i,"Variant_Classification"]) )
#      mut_lab_1 = paste( mut.dat[i,"locus"], ", ", var_info, sep="" )
      
#  tabcols = c("SC.ix", "C.ix", "total_qc", "total_qs", "SC_Aq_d", "SC_Aq_a", "C_Aq", "A1.ix", "A2.ix")

      mtext( text=mut_lab_1, side=3, line= 0.5, adj=0, cex=par("cex"), font=1 )
      mtext( text=mut_lab_2, side=3, line=-0.3, adj=0, cex=par("cex"), font=1 )
   }
   
   par(mar=dmar)
}
