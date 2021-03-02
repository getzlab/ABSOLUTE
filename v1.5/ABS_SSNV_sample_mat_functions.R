
reject_mutation_classes = function()
{
   return( c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region", "De_novo_Start_InFrame") )
}

fillin.samps = function( res_mat, gene_list, sample_list )
{
#   missing_samples = setdiff(sample_list, colnames(res_mat) )
   xx = setdiff(colnames(res_mat), sample_list )
   if( length(xx)>0 ) { stop("missing samples??") }

   mat = matrix( 0, nrow=length(gene_list), ncol=length(sample_list))
   colnames(mat) = sample_list
   rownames(mat) = gene_list
   mat[, colnames(res_mat)] = res_mat

   return(mat)
}

get_variant_class_levels = function()
{
   return( c("Synonymous", "In_frame_Indel", "Other_non_syn.", "Missense", "Splice_Site", "Frame_Shift", "Nonsense", "miRNA", "Silent", "lincRNA", "homozygous deletion", "amplification", "high-level amplification", "Missense_COSMIC_site" ) )
}

get_variant_classification_number = function(MAF)
{
   Variant_Classification = MAF[["Variant_Classification"]]
   hotspot.ix = MAF[["Number of times codon change is in COSMIC"]]

   tmp <- lapply( list(c("_Mutation", ""), 
                       c("Silent", "Synonymous"), 
                       c("Splice_Site.*", "Splice_Site"), 
                       c("Nonstop|De_novo_Start.*|Start_.*|Translation_.*|Read\\-through.*", "Other_non_syn."), 
                       c("In_frame.*", "In_frame_Indel"), 
                       c("Frame_Shift.*", "Frame_Shift"), 
                       c("3'\\-?UTR|5'\\-?UTR|3'\\-?Flank|5'\\-?Flank|IGR|Intron", "Silent")), 
                  function(x) Variant_Classification <<- sub(x[1], x[2], as.character(Variant_Classification), ignore.case=TRUE))

   Variant_Classification[ Variant_Classification == "Missense" & hotspot.ix ] = "Missense_COSMIC_site"

   Variant_Classification_Num <- as.numeric(factor(Variant_Classification, levels=get_variant_class_levels() ))

   return(Variant_Classification_Num)
}

get_mut_count_matrix = function( agg_MAF, gene_list, sample_list )
{
   final_analysis_set = agg_MAF
   final_analysis_set$Hugo_Symbol <- factor(final_analysis_set$Hugo_Symbol, levels=gene_list)
   final_analysis_set$count <- 1

## group recurrent mutations in the same gene/sample
   grouped <- aggregate(count ~ Hugo_Symbol + pair_id, sum, na.rm=TRUE, data=final_analysis_set)
   res_mat <- xtabs(count ~ Hugo_Symbol + pair_id, data=grouped, drop.unused.levels = FALSE)
   mat = fillin.samps( res_mat, gene_list, sample_list )
   return(mat)
}

aggregate_numeric_MAF_column = function( agg_MAF, gene_list, sample_list, colname, agg_func )
{
   final_analysis_set = agg_MAF
   final_analysis_set$Hugo_Symbol <- factor(final_analysis_set$Hugo_Symbol, levels=gene_list)
   final_analysis_set$count <- final_analysis_set[,colname]

## group recurrent mutations in the same gene
   grouped <- aggregate(count ~ Hugo_Symbol + pair_id, agg_func, na.rm=TRUE, data=final_analysis_set)
   res_mat <- xtabs(count ~ Hugo_Symbol + pair_id, data=grouped, drop.unused.levels = FALSE)
   mat = fillin.samps( res_mat, gene_list, sample_list )
   return(mat)
}

aggregate_boolean_MAF_column = function( agg_MAF, gene_list, sample_list, annot, agg_func )
{
   if( !(annot %in% colnames(agg_MAF) ) ) { stop("Invalid annot given") }

   final_analysis_set = agg_MAF
   final_analysis_set$Hugo_Symbol <- factor(final_analysis_set$Hugo_Symbol, levels=gene_list)
   final_analysis_set$res <- as.numeric(factor(as.integer(final_analysis_set[[annot]]), levels=c(0,1) )) 

   agg <- aggregate(res ~ Hugo_Symbol + pair_id, agg_func, na.rm=TRUE, data=final_analysis_set)
   res_mat <- xtabs(res ~ Hugo_Symbol + pair_id, data=agg, drop.unused.levels = FALSE)
   mat = fillin.samps( res_mat, gene_list, sample_list )

   return(mat)
}


make_multi_VCF = function( MAF, gene_list, ordered_sample_names, annot_names, TF_annot )
{
   multi_VCF = array(0, dim=c(length(gene_list), length(ordered_sample_names), length(annot_names)))
   dimnames(multi_VCF)[[1]] = gene_list
   dimnames(multi_VCF)[[2]] = ordered_sample_names
   dimnames(multi_VCF)[[3]] = annot_names

   if( nrow(MAF)==0 )
   {
      return( multi_VCF )
   }

   multi_VCF[,,"var_class"] = aggregate_numeric_MAF_column( MAF, gene_list, ordered_sample_names, "Variant_Classification_Num", max )

#   multi_VCF[,,"ssnv.biallelic"] = get_mut_count_matrix( MAF, gene_list, ordered_sample_names )

   multi_VCF[,,"clinical actionability priority"] = aggregate_numeric_MAF_column( MAF, gene_list, ordered_sample_names, "clinical actionability priority", min )

   if( "underpowered_in_paired_sample" %in% colnames(MAF) )
   {
      nix = is.na(MAF[,"underpowered_in_paired_sample"])
      MAF[nix,"underpowered_in_paired_sample"] = 0    ## assume SCNAs are powered
      multi_VCF[,,"underpowered_in_paired_sample"] = aggregate_numeric_MAF_column( MAF, gene_list, ordered_sample_names, "underpowered_in_paired_sample", max)
   }

   for( i in 1:length(TF_annot)) 
   {
      multi_VCF[ ,, TF_annot[i] ] = aggregate_boolean_MAF_column( MAF, gene_list, ordered_sample_names, TF_annot[i], max )
   }

   return(multi_VCF)
}


## create a set of multidimensional VCFs from a MAF - dims give annotations for obs muts
get_multi_VCF = function( agg_MAF, gene_list, ordered_sample_names )
{
   required_cols = c("pair_id", "Hugo_Symbol" )
   if( !all( required_cols %in% colnames(agg_MAF) ) ) { stop("Missing required agg_MAF cols!") }

   TF_annot = c("homozygous.ix", "ssnv.ge2.ix"  )
#   annot_names = c("var_class", "COSMIC_site_count", "ssnv.biallelic", "clinical actionability priority", TF_annot )
   annot_names = c("var_class", "clinical actionability priority", TF_annot )
   if( "underpowered_in_paired_sample" %in% colnames(agg_MAF) )
   {
      annot_names = c(annot_names, "underpowered_in_paired_sample" ) 
   }

## fill in missing optional cols
   for( i in 1:length(TF_annot) )
   {
      if( !(TF_annot[i] %in% colnames(agg_MAF)) )
      {
        ## Fill in missing cols with FALSE
         agg_MAF = cbind( agg_MAF, rep(FALSE, nrow(agg_MAF)) )
         names(agg_MAF)[ncol(agg_MAF)] = TF_annot[i]
      }

     ## Set NAs to FALSE
      nix = is.na( agg_MAF[, TF_annot[i] ] )
      agg_MAF[ nix, TF_annot[i] ] = FALSE
   }


   filtered_MAF = agg_MAF[ agg_MAF[,"Hugo_Symbol"] %in% gene_list,]
   nix = is.na(filtered_MAF[,"Number of times codon change is in COSMIC"])
   filtered_MAF[nix,"Number of times codon change is in COSMIC"] = 0
   filtered_MAF[["Variant_Classification_Num"]] = get_variant_classification_number( filtered_MAF )

  ## reorder so most 'severe' mutations are 1st
   o.ix = order( filtered_MAF[, c("pair_id", "Hugo_Symbol", "Variant_Classification_Num")], decreasing=TRUE )
   filtered_MAF = filtered_MAF[ o.ix,] 

## find multiple muts in same gene/sample
   second.mut.ix = duplicated( filtered_MAF[, c("pair_id", "Hugo_Symbol")] )

# more severe mut will be primary..
   primary_MAF = filtered_MAF[ !second.mut.ix,] 
   secondary_MAF = filtered_MAF[ second.mut.ix, ]
#   third.mut.ix = duplicated( second_filtered_MAF[, c("pair_id", "Hugo_Symbol")] )

   primary_multi_VCF = make_multi_VCF( primary_MAF, gene_list, ordered_sample_names, annot_names, TF_annot )
   secondary_multi_VCF = make_multi_VCF( secondary_MAF, gene_list, ordered_sample_names, annot_names, TF_annot )

   return( list("primary_multi_VCF"=primary_multi_VCF, "secondary_multi_VCF"=secondary_multi_VCF) )
}


## Co-mut style plot with muts colored by variant type (indel, nonsense, etc)
plot_mutmat = function( gene_matrix, tracks_mat, track_colors, mut_rate_color, mut_rates, max_mut_rate, label, plot_legend=FALSE )
{
   genematrix = gene_matrix[,,1]
   hom_ssnv_matrix = gene_matrix[,,2]
   mut_count_mat = gene_matrix[,,3]

   mutation_rates=TRUE
   plot_tracks=TRUE
  
   N_genes = nrow(genematrix)
   N_samps = ncol(genematrix)

   par(mar=c(0.5,7,0,2), las=1)

   mut_type_colors =  col=c("grey94", brewer.pal(7, "Set1")[c(3,6,7,2,4,5,1)])
   image( t(genematrix), col=mut_type_colors, zlim=c(0,7), xlab="", ylab="", axes=FALSE)

   yc = par("usr")[c(3,4)]
   yrng = yc[2] - yc[1]
   axis( side=2, at= yc[1] + yrng/N_genes * c(1:N_genes-0.5), labels=rownames(genematrix), font=3, las=1,
         cex.axis=1.0, tck=-0.005, mgp=c(3, 0.3, 0) )

   xc = par("usr")[c(1,2)]
   xrng = xc[2] - xc[1]

## Homozygous mutation
   h.ix = which(hom_ssnv_matrix==2, arr.ind=T)
   points( xc[1] + xrng/N_samps * (h.ix[,2]-0.5), 
           yc[1] + yrng/N_genes * (h.ix[,1]-0.5), col="white", pch=16, cex=0.6 )

## > 1 mutation in gene/samp
   m.ix = which(mut_count_mat>1, arr.ind=T)
   points( xc[1] + xrng/N_samps * (m.ix[,2]-0.5), 
           yc[1] + yrng/N_genes * (m.ix[,1]-0.5), col="white", pch="+", cex=0.6 )

## sample categories
   if( plot_tracks ) {
      image( tracks_mat, col=track_colors, xlab="", ylab="", axes=FALSE)
   }


   if( mutation_rates )
   {
	## Mutation rates
	par(mar=c(.25, 7, 1, 2), las=1)
#	if (verbose) cat("  mutations rates\n")
	plot(1, ylim=c(0,max_mut_rate), type="n", axes=F, xlab="", ylab="")
	par(usr=c(0, ncol(genematrix), par("usr")[c(3,4)]), lwd=ifelse(ncol(genematrix)>200, .2, .5))

	barplot(mut_rates, col=mut_rate_color, axes=FALSE, add=TRUE, names.arg=rep("", length(mut_rates)), border="grey94", space=0)
	par(lwd=1)
	axis(2, cex.axis=1, line=.3)
#	if (any(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max_mut_rate)) {
 # 	   text(which(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max_mut_rate)-.5, max_mut_rate*.9, cex=.8, col="white", labels=round(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6))[which(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max_mut_rate)], srt=90)
#        }
	mtext("# mutations/Mb", side=2, cex=.7, line=2.7, las=0)

	mtext(label, side=3, cex=0.7, line=-0.3, font=2, adj=0 )

   }

   if( plot_legend ) 
   {
	plot(0, type="n", axes=FALSE, xlab="", ylab="")
	par(xpd=NA)
	smartlegend("right", "top", c("Syn.", "Missense", "Splice site", "Nonsense", "Frame shift", "In frame indel", "Other non syn."),
			fill=brewer.pal(7, "Set1")[c(3,2,4,1,5,6,7)], cex=.9, ncol=2, bty = "n", inset=0)
	par(xpd=FALSE)
   }
}


plot_mut_rates_and_sample_tracks = function( tracks_mat, tracks_colors, mut_colors, mut_rates )
{
   N_samps = ncol(mut_rates)
   rmar=0.2
   par(mar=c(0.5,7,0,rmar), las=1)
   image( tracks_mat, col=track_colors, xlab="", ylab="", axes=FALSE)

## barplot
   par(mar=c(.25, 7, 2, rmar), las=1)
   plot(1, ylim=c(0,max_mut_rate), type="n", axes=F, xlab="", ylab="")
   par(usr=c(0, N_samps, par("usr")[c(3,4)]), lwd=ifelse(N_samps>200, .2, .5))

   barplot(mut_rates, col=mut_colors, axes=FALSE, add=TRUE, names.arg=rep("", ncol(mut_rates)), border=NA, space=0)
   par(lwd=1)
   axis(2, cex.axis=1, line=.3)
   mtext("SSNVs / Mb", side=2, cex=.7, line=2.7, las=0)
}


# works from output .maf files..
# much slower than outputting from segobj.list.
aggregate_ABS_MAF_files = function( MAF_FNs )
{
   combined_MAF = data.frame()
   MAF_list = list()
   for( i in 1:length(MAF_FNs))
   {
      MAF = read.delim(MAF_FNs[i], row.names = NULL, stringsAsFactors = FALSE, 
                        check.names = FALSE, na.strings = c("NA", "---"),
                        blank.lines.skip=TRUE, comment.char="#")

      MAF_list[[i]] = MAF
   }

   return(MAF_list)
}


get_MAF_list_from_called_seglist_obj = function( called_segobj_list )
{
 MAF_list = vector( length=length(called_segobj_list), mode="list")
 nix = rep( FALSE, length(called_segobj_list))

 for (i in seq_along(called_segobj_list)) 
 {
    called_segobj = called_segobj_list[[i]] 
    names(MAF_list)[i] = called_segobj[["sample.name"]]
    
    ## MAF
    if (!is.null(called_segobj$mode.res$modeled.muts) &&   ## need the short-circuit
         nrow(called_segobj$mut.cn.dat) > 0 ) 
    {
      modeled = called_segobj$mode.res$modeled.muts[[1]]
      mut_dat = cbind(called_segobj$mut.cn.dat, modeled)      

      SSNV_ccf_dens = called_segobj[["mode.res"]][["SSNV.ccf.dens"]][1,,]
      MAF_list[[i]] = data.frame( "pair_id"=names(called_segobj_list)[i], mut_dat, SSNV_ccf_dens, check.names=FALSE, stringsAsFactors=FALSE)
    }
    else{ nix[i] = TRUE }
  }

  if( any(nix) ) { MAF_list = MAF_list[!nix] }

  return(MAF_list)
}

aggregate_sample_MAF_list = function( MAF_list )
{
   cols = colnames(MAF_list[[1]] )

   if( length(MAF_list) > 1 )
   {
      for( i in 2:length(MAF_list))
      {
         cols = intersect(cols, colnames(MAF_list[[i]]))
      }
   }

   combined_MAF = data.frame()
   for( i in 1:length(MAF_list))
   {
      combined_MAF = rbind( MAF_list[[i]][,cols], combined_MAF ) 
   }

   return(combined_MAF)
}



