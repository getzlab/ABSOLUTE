get_sample_background_rate_convolved_event_significance = function( event_mat, sample_bg_rate, PDF_OUT_FN, TSV_OUT_FN )
{
   sbg = sample_bg_rate
   
   ugenes = rownames(event_mat)
   cols = c( "P_val", "Q_val", "N_alt_class", "N_alt_other")
   gmat = matrix( 0, nrow=length(ugenes), ncol=length(cols))
   rownames(gmat)=ugenes
   colnames(gmat)=cols
   gmat[, "N_alt_class"] = rowSums( event_mat==2 )
   gmat[, "N_alt_other"] = rowSums( event_mat==1 )

   for( i in 1:nrow(gmat) )
   {
      if( gmat[i,"N_alt_class"] == 0 ) { gmat[i,"P_val"] = 1; next }  ## if the event was never observed
   
      s.ix = which(event_mat[i,] > 0)
      ev_samp_bgs = sbg[s.ix]
      comb_pr = c(1-ev_samp_bgs[1], ev_samp_bgs[1])             ## if the event was observed once
   
      if( length(ev_samp_bgs) > 1 )                             ## if the event was observed > once
      {
         for( j in 2:length(ev_samp_bgs) )
         {
            vp = c( 1-ev_samp_bgs[j], ev_samp_bgs[j])
            comb_pr = convolve(comb_pr, rev(vp), type = "o")
         }
      }
      
     ## P-val is Pr(getting >= obs hom_muts (k) under null comb_pr)
      k = gmat[i,"N_alt_class"]
      gmat[i,"P_val"] = sum( comb_pr[ c(k:(length(comb_pr)-1)) + 1 ] )
   }
   
   gmat[,"Q_val"] = gmat[,"P_val"] * nrow(gmat)  ## simple Bonferroni
   gmat[gmat[,"Q_val"] > 1, "Q_val"] = 1
   
   return(gmat)
}  




homozygous_mutsig = function( agg_MAF, MAF_list, called.segobj.list, PDF_OUT_FN, TSV_OUT_FN, ordered_sample_names=NA )  
{
   samp_frac_LOH = rep(NA, length(called.segobj.list) )
   names(samp_frac_LOH) = names(called.segobj.list) 

   for( i in 1:length(called.segobj.list))
   {
      abs_segs = get_allelic_abs_CN_seg_dat_from_allelic_CAPSEG_obj( called.segobj.list[[i]] )
      samp_frac_LOH[i] = sum(abs_segs[,"LOH"] * abs_segs[,"W"])
   }

   samp_SSNV_mat = matrix( NA, nrow=length(MAF_list), ncol=3 )
   rownames(samp_SSNV_mat) = names(MAF_list) 
   hom_class_cols = c("wt0.ix", "subclonal_wt0.ix")
   for( i in 1:length(MAF_list) )
   {
      h.ix = apply( MAF_list[[i]][, hom_class_cols ], 1, any, na.rm=TRUE )
      n.ix =  is.na(h.ix)
      h.ix[n.ix] = FALSE
      samp_SSNV_mat[i, ] = c(sum(h.ix), sum(!h.ix), sum(h.ix) / length(h.ix) )
   }

   agg_MAF = agg_MAF[ !(agg_MAF[, "Variant_Classification"] %in% reject_mutation_classes()), ]

   if( is.na(ordered_sample_names))
   {
      ordered_sample_names = sort( unique( agg_MAF[,"pair_id"] ))
   }
   gene_list = sort( unique( agg_MAF[,"Hugo_Symbol"] ))
   mut_mats = get_multi_VCF( agg_MAF, gene_list, ordered_sample_names )
   samp_SSNV_mat = samp_SSNV_mat[ordered_sample_names,]
   samp_frac_LOH = samp_frac_LOH[ordered_sample_names]
   samp_frac_SSNV_hom = samp_SSNV_mat[,3]

   homozygous_SSNV_mat = mut_mats[,,2]  ## values of 2 indicate hom mut

   gmat_LOH = get_sample_background_rate_convolved_event_significance( homozygous_SSNV_mat, samp_frac_LOH, PDF_OUT_FN, TSV_OUT_FN )

   gmat_frac = get_sample_background_rate_convolved_event_significance( homozygous_SSNV_mat, samp_frac_SSNV_hom, PDF_OUT_FN, TSV_OUT_FN )


   pdf(PDF_OUT_FN, 10, 10 )
   par(mfrow=c(2,2))
   par(bty="n", las=1)

   plot( samp_frac_LOH, samp_frac_SSNV_hom, main="", xlab="Fraction of genome at LOH / sample", ylab="Fraction of homozygous SSNVs / sample", pch=16, cex=0.5 )

   samp_frac_SSNV_hom_95 = cbind(  qbeta( 0.025, samp_SSNV_mat[,1] + 1, samp_SSNV_mat[,2] + 1), 
                                   qbeta( 0.975, samp_SSNV_mat[,1] + 1, samp_SSNV_mat[,2] + 1) )

   segments( y0=samp_frac_SSNV_hom_95[,1], y1=samp_frac_SSNV_hom_95[,2], x0=samp_frac_LOH, x1=samp_frac_LOH  ) 

   frame()
   qq_pval( gmat_LOH[,"P_val"], qvalues=gmat_LOH[,"Q_val"], genes=rownames(gmat_LOH) )
   qq_pval( gmat_frac[,"P_val"], qvalues=gmat_frac[,"Q_val"], genes=rownames(gmat_frac) )
   
   dev.off()


   gmat = gmat_LOH
   o.ix = order(gmat[,"P_val"])
   gmat = gmat[o.ix,]

   write.table( gmat, file=TSV_OUT_FN, quote=FALSE, sep="\t")

}


ABS_mutsig_calc_and_plot = function( mutclass_mat, nonsil_mutclass_mat, samp_frac_genome_at_risk, PDF_OUT_FN, TSV_OUT_FN, XLAB, YLAB )
{
   samp_SSNV_mat = cbind( colSums(mutclass_mat>1), colSums(mutclass_mat==1))
   samp_SSNV_mat = cbind( samp_SSNV_mat, samp_SSNV_mat[,1] / rowSums(samp_SSNV_mat) )
   samp_frac_SSNV_class = samp_SSNV_mat[,3]

   gmat_f_genome_at_risk = get_sample_background_rate_convolved_event_significance( nonsil_mutclass_mat, samp_frac_genome_at_risk, PDF_OUT_FN, TSV_OUT_FN )

   gmat_f_SSNVs = get_sample_background_rate_convolved_event_significance( nonsil_mutclass_mat, samp_frac_SSNV_class, PDF_OUT_FN, TSV_OUT_FN )


   pdf(PDF_OUT_FN, 10, 10 )
   par(mfrow=c(2,2))
   par(bty="n", las=1)

   plot( samp_frac_genome_at_risk, samp_frac_SSNV_class, main="", xlab=XLAB, ylab=YLAB, pch=16, cex=0.5 )

   samp_frac_SSNV_class_95 = cbind(  qbeta( 0.025, samp_SSNV_mat[,1] + 1, samp_SSNV_mat[,2] + 1), 
                                   qbeta( 0.975, samp_SSNV_mat[,1] + 1, samp_SSNV_mat[,2] + 1) )

   segments( y0=samp_frac_SSNV_class_95[,1], y1=samp_frac_SSNV_class_95[,2], x0=samp_frac_genome_at_risk, x1=samp_frac_genome_at_risk  ) 

   frame()
   qq_pval( gmat_f_genome_at_risk[,"P_val"], qvalues=gmat_f_genome_at_risk[,"Q_val"], genes=rownames(gmat_f_genome_at_risk) )
   qq_pval( gmat_f_SSNVs[,"P_val"], qvalues=gmat_f_SSNVs[,"Q_val"], genes=rownames(gmat_f_SSNVs) )
   
   dev.off()


   gmat = gmat_f_SSNVs
   o.ix = order(gmat[,"P_val"])
   gmat = gmat[o.ix,]

   write.table( gmat, file=TSV_OUT_FN, quote=FALSE, sep="\t")
}



ABS_class_mutsig = function( agg_MAF, called.segobj.list, out.base )
{
   ordered_sample_names = sort( unique( agg_MAF[,"pair_id"] ))
   gene_list = sort( unique( agg_MAF[,"Hugo_Symbol"] ))
   mut_mats = get_multi_VCF( agg_MAF, gene_list, ordered_sample_names )

   homozygous_SSNV_mat = mut_mats[,,2]  ## values of 2 indicate hom mut
   ge2_SSNV_mat = mut_mats[,,4]  
#   double_hit = mut_mats[,,3]
#   homozygous_SSNV_mat[ double_hit > 1 ] = 2

## Now remove synon / non-coding muts and regenerate...
   agg_MAF = agg_MAF[ !(agg_MAF[, "Variant_Classification"] %in% reject_mutation_classes()), ]
   gene_list = sort( unique( agg_MAF[,"Hugo_Symbol"] ))
   nonsil_mut_mats = get_multi_VCF( agg_MAF, gene_list, ordered_sample_names )

   nonsil_homozygous_SSNV_mat = nonsil_mut_mats[,,2]
   nonsil_ge2_SSNV_mat = nonsil_mut_mats[,,4]
#   nonsil_double_hit = nonsil_mut_mats[,,3]
# count >1 hit genes as homozygous for significance calc
# Actually this is a very bad idea - get flooded with high mutation rate genes  CSMD3, etc
#   nonsil_homozygous_SSNV_mat[ nonsil_double_hit > 1 ] = 2

   samp_frac_gain = rep(NA, length(called.segobj.list) )
   samp_frac_LOH = rep(NA, length(called.segobj.list) )
   names(samp_frac_gain) = names(called.segobj.list) 
   names(samp_frac_LOH) = names(called.segobj.list) 

   for( i in 1:length(called.segobj.list))
   {
      abs_segs = get_allelic_abs_CN_seg_dat_from_allelic_CAPSEG_obj( called.segobj.list[[i]] )
      samp_frac_gain[i] = sum((abs_segs[,"modal.a1"]>1 | abs_segs[,"modal.a2"]>1) * abs_segs[,"W"])
      samp_frac_LOH[i] = sum(abs_segs[,"LOH"] * abs_segs[,"W"])
   }
   samp_frac_gain = samp_frac_gain[ordered_sample_names]
   samp_frac_LOH = samp_frac_LOH[ordered_sample_names]


   ABS_mutsig_calc_and_plot( homozygous_SSNV_mat, nonsil_homozygous_SSNV_mat, samp_frac_LOH, PDF_OUT_FN=paste(out.base, "homozygous_mutsig.pdf", sep=""), TSV_OUT_FN=paste( out.base, "homozygous_mutsig.txt", sep=""), XLAB="Fraction of genome at LOH / sample", YLAB="Fraction of homozygous SSNVs / sample" )
  
   ABS_mutsig_calc_and_plot( ge2_SSNV_mat, nonsil_ge2_SSNV_mat, samp_frac_gain, PDF_OUT_FN=paste( out.base, "ge2_mutsig.pdf",sep=""), TSV_OUT_FN=paste( out.base, "ge2_mutsig.txt", sep=""), XLAB="Fraction of genome at > 1 AS copy / sample" , YLAB="Fraction of ge2 SSNVs / sample" )

}

