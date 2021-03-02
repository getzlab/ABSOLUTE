build_gene_GR_data = function()
{
   #genes_FN = "refGene.hg19.20100825.sorted.txt"
   #genes_FN = "/xchip/cga/reference/annotation/db/ucsc/hg19/gene_table.txt"
#   genes_FN = "/xchip/gistic/variables/hg19/rg_20120227_dump.txt"
#  gene_dat = read.delim(genes_FN, check.names=FALSE, stringsAsFactors=FALSE )

   #data(refgene)
   load("/xchip/tcga/Tools/absolute/releases/v1.5/data/refgene.hg19.genes.RData")   
   gene_dat = refgene 

   gene_dat[,"chr"] =  gsub("chr", "", gene_dat[,"chr"] )
   gene_footprints = GRanges( gene_dat[,"chr"], IRanges(gene_dat[,"start"], gene_dat[,"end"]) )
#   colnames(gene_dat)[3] = "name"   ## symb -> name
   names(gene_footprints) = gene_dat[,3]
   
   #name	longname	chr	strand	tx_start	tx_end
   #gene_footprints = GRanges( gene_dat[,"chr"], IRanges(gene_dat[,"tx_start"], gene_dat[,"tx_end"]) )
   #gene_footprints = GRanges( gene_dat[,"chr"], IRanges(gene_dat[,"gene_start"], gene_dat[,"gene_end"]) )

  return( gene_footprints )
}

load_GISTIC_peak_GR_data = function( GISTIC_peaks_fn, reg_tab=NA )
{ 
   GISTIC_peaks = read.delim( GISTIC_peaks_fn, check.names=FALSE, stringsAsFactors=FALSE, header=0, row.names=1, nrows=5 )
   GISTIC_peaks = GISTIC_peaks[,-ncol(GISTIC_peaks)]   ## extra col
   peak_strs = as.character( GISTIC_peaks["wide peak boundaries",] )
   
   peak_mat =  matrix( unlist(strsplit( peak_strs, ":|-")), ncol=3, byrow=TRUE )
   peak_mat[,1] = gsub( "chr", "", peak_mat[,1] )
   
   peak_GR = GRanges( peak_mat[,1], IRanges(as.integer(peak_mat[,2]), as.integer(peak_mat[,3])) )
   names(peak_GR) =  GISTIC_peaks["wide peak boundaries",]

# try to rename regs
   if( !is.na(reg_tab))
   {
      r.ix = match( names(peak_GR), reg_tab[,"Peak region"] )
      if( any(is.na(r.ix)) ) { stop("unmatched reg name!") }

      reg_names = reg_tab[r.ix, "Peak Name"]
      names(peak_GR) = reg_names
   }

   return(peak_GR)
} 


compute_focality_score = function( segtab, mode )
{
   if( !(mode %in% c("amp", "del") ) ) { stop() }

   cn = segtab[,"rescaled_total_cn"]
   WD = 0.1
   br = seq( min(cn), max(cn)+WD, by=WD)

   bin = rep(NA, length(cn)) 
   for (i in 1:length(cn) )
   {
      bin[i] <- max(which(br <= cn[i]))
   }

   levels = sort(unique(bin))

   if( mode == "del" ) { levels = rev(levels) }

   genome_frac = rep(NA, length(levels)) 
   for( i in 1:length(levels) )
   {
      genome_frac[i] = sum(segtab[bin==levels[i],"W"])
   }  

   cum_genome_frac = cumsum(c(0,genome_frac) )
   focality = rep(NA, nrow(segtab))
   for( i in 1:length(levels) )
   {
     focality[bin==levels[i]] = cum_genome_frac[i] 
   }  

   return(focality)
}


## Calls all segments
call_genome_wide_ABSOLUTE_SCNAs = function( ABS.dat )
{
   segtab = AllelicGetAbsSegDat(ABS.dat)

   seg_amp_focality = compute_focality_score(segtab, "amp" )
   seg_del_focality = compute_focality_score(segtab, "del" )

   total.only.ix = is.na( segtab[,"HZ"] )
   hzdel.ix = rep( FALSE, nrow(segtab) )


#   hzdel.ix[!total.only.ix] =  (segtab[!total.only.ix,"total_HZ"] | segtab[!total.only.ix,"HZ"] | segtab[!total.only.ix, "SC_HZ" ]) & segtab[!total.only.ix,"corrected_total_cn"] < 0.75 & seg_del_focality[!total.only.ix] > 0.99
#   hzdel.ix[!total.only.ix] =  (segtab[!total.only.ix,"total_HZ"] | segtab[!total.only.ix,"HZ"] | segtab[!total.only.ix, "SC_HZ" ]) & segtab[!total.only.ix,"corrected_total_cn"] < 0.75 & seg_del_focality[!total.only.ix] > 0.99
#   hzdel.ix[!total.only.ix] = segtab[!total.only.ix,"rescaled_total_cn"] < 0.25 & seg_del_focality[!total.only.ix] > 0.995
#   hzdel.ix[ total.only.ix] =  segtab[ total.only.ix,"total_HZ"]
#   hzdel.ix[ total.only.ix] = segtab[total.only.ix,"rescaled_total_cn"] < 0.25 & seg_del_focality[total.only.ix] > 0.995

## This works OK - hom dels in CDKN2A can be up to 11MB! (0.5% genome)
   hzdel.ix = segtab[,"rescaled_total_cn"] < 0.25 & seg_del_focality > 0.995

   del = as.logical(hzdel.ix)

#print(files[i])
#if(  length(grep("PB0271", files[i])) > 0 ) { print(segtab[49:54,]); stop() }

#if( any( seg_del_focality[!total.only.ix] > 0.99 & segtab[!total.only.ix,"corrected_total_cn"] < 0.75 & !hzdel.ix[!total.only.ix] ) ) { stop() }


## AMPS
   ploidy = ABS.dat[["mode.res"]][["mode.tab"]][1,"genome mass"]
   purity = ABS.dat[["mode.res"]][["mode.tab"]][1,"alpha"]
   WGD = as.integer(ABS.dat[["mode.res"]][["mode.tab"]][1,"WGD"])

#   amp = segtab[,"corrected_total_cn"] >= 7.0
   amp = seg_amp_focality > 0.98 & segtab[,"rescaled_total_cn"] >= 5.0

   called_segtab = cbind( segtab, "amp.call"=amp, "del.call"=del, "amp_foc"=seg_amp_focality, "del_foc"=seg_del_focality, "ploidy"=ploidy, "purity"=purity, "WGD"=WGD )

   return( called_segtab )

}

## Convert calls on segs to calls on defined regions e.g. peaks, genes
call_regions_ABSOLUTE_SCNAs = function( segtab, regions )
{
   seg_GRs = GRanges( segtab[,"Chromosome"], IRanges(segtab[,"Start.bp"], segtab[,"End.bp"]) )

# gr.findoverlaps(gr1, gr2) is a faster and less clumsy reimplementation of GRanges findOverlaps
#   reg_to_seg = gr.match( regions, seg_GRs )  ## only returns 1st match!!  not cool
   reg_to_seg = gr.findoverlaps( regions, seg_GRs )

## Need to resolve multiple matches from regions to genomic segs.   Break ties using CN (highest match for amps, lowest for dels)
   amp.seg.ix = rep(NA, length(regions) )
   del.seg.ix = rep(NA, length(regions) )
   for( i in 1:length(regions) )
   {
      ix = which( reg_to_seg$query.id == i )
      if( length(ix) == 0 ) { next }
      CN_vals = segtab[ reg_to_seg$subject.id[ ix ], "corrected_total_cn"]
      amp.seg.ix[i] = reg_to_seg$subject.id[ ix[ which.max(CN_vals) ] ]
      del.seg.ix[i] = reg_to_seg$subject.id[ ix[ which.min(CN_vals) ] ]
   }
   
# Some regions may not hit any segs because the reg is between two segs
   missing.ix = which( is.na(amp.seg.ix) )
   for( i in missing.ix )
   {
      dd = distance( regions[i], seg_GRs )
# indices of 2 closest seg matches
      ix2 = order(dd, decreasing=FALSE, na.last=TRUE)[c(1,2)]
 
      CN_vals = segtab[ ix2, "corrected_total_cn"]
      amp.seg.ix[i] = ix2[ which.max(CN_vals) ] 
      del.seg.ix[i] = ix2[ which.min(CN_vals) ] 
   }

# create segtabs containing unique 'best' segment matching each input region
   amp.reg.segtab = segtab[ amp.seg.ix, ]
   del.reg.segtab = segtab[ del.seg.ix, ]

   return( list("amp.reg.segtab"=amp.reg.segtab, "del.reg.segtab"=del.reg.segtab) )
}

genotype_transcript_SCNAs_in_called_ABS_files = function( ABS_BASE_DIR, regs, SCNA_thresholds, analyst_id="SLC" )
{
   fn_exts = paste(".ABSOLUTE.", analyst_id, ".called.RData", sep="")
   files = grep(  fn_exts, dir(ABS_BASE_DIR, full.names=FALSE), value=TRUE)
   snames = gsub( fn_exts, "", files )

   amp_ev_mat = matrix( NA, nrow=length(regs), ncol=length(snames) )
   del_ev_mat = matrix( NA, nrow=length(regs), ncol=length(snames) )

   rownames(amp_ev_mat) = names(regs)
   rownames(del_ev_mat) = names(regs)
   colnames(amp_ev_mat) = colnames(del_ev_mat) = snames

### Debugging:
#   i = 24
#   i = grep( "PB0274-TM", files )
#   load( file.path(ABS_BASE_DIR, files[i]) ) 
#   ABS.dat = seg.obj
#   called.segtab = call_genome_wide_ABSOLUTE_SCNAs( ABS.dat )
#   del.reg.segtab = call_regions_ABSOLUTE_SCNAs( called.segtab, regs ) [["del.reg.segtab"]]

#   called.segtabs = list()
#   amp.segtabs = list()
#   del.segtabs = list()
#   for( i in 1:length(files) )
   res = foreach( i = 1:length(files)) %dopar%
   {
      load( file.path(ABS_BASE_DIR, files[i]) ) 
      ABS.dat = seg.obj

      called.segtab = call_genome_wide_ABSOLUTE_SCNAs( ABS.dat )
      res = call_regions_ABSOLUTE_SCNAs( called.segtab, regs ) 

      del.reg.segtab = res[["del.reg.segtab"]]
      amp.reg.segtab = res[["amp.reg.segtab"]]

      cat(".")
      return( list("called.segtab"=called.segtab, "amp.segtab"=amp.reg.segtab, "del.segtab"=del.reg.segtab) )

#      called.segtabs[[i]] = called.segtab
#      amp.segtabs[[i]] = amp.reg.segtab
#      del.segtabs[[i]] = del.reg.segtab
   }

   called.segtabs = lapply( res, "[[", "called.segtab" )
   amp.segtabs = lapply( res, "[[", "amp.segtab" )
   del.segtabs = lapply( res, "[[", "del.segtab" )

   for( i in 1:length(files) )
   {
      amp_ev_mat[,i] = amp.segtabs[[i]][,"amp.call"]
      del_ev_mat[,i] = del.segtabs[[i]][,"del.call"]
   }
   cat("done\n")

   names(called.segtabs) = names(amp.segtabs) = names(del.segtabs) = snames
   SCNA_event_dat = list("amp_ev_mat"=amp_ev_mat, "del_ev_mat"=del_ev_mat, "amp_regs"=regs, "del_regs"=regs )  

# build 3D array: genes X samples X annotations for amp SCNAs  X  refgene transcripts (genes)
   genes = names(regs)  ## exactly the same as del_regs
   N_samps = length(called.segtabs)
   amp.gene.data = array( NA, dim=c(length(genes), N_samps, 5) )
   dimnames(amp.gene.data)[[1]] = genes
   dimnames(amp.gene.data)[[2]] = names(called.segtabs)
   cols = c("amp_foc", "corrected_total_cn", "rescaled_total_cn", "amp.call", "length")
   dimnames(amp.gene.data)[[3]] = cols
   
   for( i in 1:length(cols) )
   {
      amp.gene.data[,,i] = matrix( unlist( lapply( amp.segtabs, "[", cols[i])), nrow=length(genes), ncol=N_samps, byrow=FALSE )
   }


# 3D array for deletions X refgene
   N_samps = length(called.segtabs)
   del.gene.data = array( NA, dim=c(length(genes), N_samps, 5) )
   dimnames(del.gene.data)[[1]] = genes
   dimnames(del.gene.data)[[2]] = names(called.segtabs)
   cols = c("del_foc", "corrected_total_cn", "rescaled_total_cn", "del.call", "length")
   dimnames(del.gene.data)[[3]] = cols
   for( i in 1:length(cols) )
   {
      del.gene.data[,,i] = matrix( unlist( lapply( del.segtabs, "[", cols[i])), nrow=length(genes), ncol=N_samps, byrow=FALSE )
   }

   SCNA_event_dat[["amp.gene.data"]] = amp.gene.data
   SCNA_event_dat[["del.gene.data"]] = del.gene.data
   
   return( list( "called.segtabs"=called.segtabs, "amp.segtabs"=amp.segtabs, "del.segtabs"=del.segtabs, "SCNA_event_dat"=SCNA_event_dat ))
}

select_samples_from_gene_SCNA_calls = function( gene_SCNA_calls, samples )
{
#[1] "called.segtabs" "amp.segtabs"    "del.segtabs"    "SCNA_event_dat"
   gene_SCNA_calls[["called.segtabs"]] = gene_SCNA_calls[["called.segtabs"]][samples]
   gene_SCNA_calls[["amp.segtabs"]] = gene_SCNA_calls[["amp.segtabs"]][samples]
   gene_SCNA_calls[["del.segtabs"]] = gene_SCNA_calls[["del.segtabs"]][samples]

   gene_SCNA_calls[["amp_ev_mat"]] = gene_SCNA_calls[["amp_ev_mat"]][,samples, drop=FALSE ]
   gene_SCNA_calls[["del_ev_mat"]] = gene_SCNA_calls[["del_ev_mat"]][,samples, drop=FALSE ]

   gene_SCNA_calls[["amp.gene.data"]] = gene_SCNA_calls[["amp.gene.data"]][ , samples, , drop=FALSE ]
   gene_SCNA_calls[["del.gene.data"]] = gene_SCNA_calls[["del.gene.data"]][ , samples, , drop=FALSE ]

   return(gene_SCNA_calls)
}




seg_focality_plots = function( PP_CALLS, SCNA_calls, drivers )
{
  ## has the genomic scnas 
   called.segtabs = SCNA_calls[["called.segtabs"]]

   wgd = PP_CALLS[ names(called.segtabs), "Genome doublings" ]
   ploidy = PP_CALLS[ names(called.segtabs), "ploidy"]
   names(ploidy) = rownames(PP_CALLS)

   amp_foc_0 =  unlist(lapply( called.segtabs[wgd==0], "[", "amp_foc" ))
   amp_foc_1 =  unlist(lapply( called.segtabs[wgd>0], "[", "amp_foc" ))
   cn_0 = unlist(lapply( called.segtabs[wgd==0], "[", "corrected_total_cn" ))
   cn_1 = unlist(lapply( called.segtabs[wgd>0], "[", "corrected_total_cn" ))

   rcn_0 = unlist(lapply( called.segtabs[wgd==0], "[", "rescaled_total_cn" ))
   rcn_1 = unlist(lapply( called.segtabs[wgd>0], "[", "rescaled_total_cn" ))

   pdf("seg_cn_vs_foc.pdf", 12, 12 )
   par(mfrow=c(2,2))
   par(  bty="n", las=1  )

# amps
   called_segs_0 = unlist(lapply(called.segtabs[wgd==0], "[", "amp.call"))
   called_segs_1 = unlist(lapply(called.segtabs[wgd>0], "[", "amp.call"))

   seg.colors_0 = rep("black", length(called_segs_0) )
   seg.colors_0[called_segs_0] = "red" 

   seg.colors_1 = rep("black", length(called_segs_1) )
   seg.colors_1[called_segs_1] = "red" 

# corrected_cn (has weighted comp mixture prior)
   plot( log(cn_0,2), amp_foc_0, main="", xlim=c(0, max(log(cn_0,2))), xlab="Log2 corrected CN", ylab="Focality", pch=".", col=seg.colors_0  )
   abline( v=log(7,2), lty=3, col=2 )

   plot( log(cn_1,2), amp_foc_1, main="", xlim=c(0, max(log(cn_1,2))), xlab="Log2 corrected CN", ylab="Focality", pch=".", col=seg.colors_1  )
   abline( v=log(7,2), lty=3, col=2 )


## rescaled_cn (linear rescale of copy-ratio using purity / ploidy)
   plot( log(rcn_0,2), amp_foc_0, main="", xlim=c(0, max(log(rcn_0,2),na.rm=TRUE)), xlab="Log2 rescaled CN", ylab="Focality", pch=".", col=seg.colors_0  )
   abline( v=log(7,2), lty=3, col=2 )

   plot( log(rcn_1,2), amp_foc_1, main="", xlim=c(0, max(log(rcn_1,2),na.rm=TRUE)), xlab="Log2 rescaled CN", ylab="Focality", pch=".", col=seg.colors_1  )
   abline( v=log(7,2), lty=3, col=2 )


## deletions
   called_segs_0 = unlist(lapply(called.segtabs[wgd==0], "[", "del.call"))
   called_segs_1 = unlist(lapply(called.segtabs[wgd>0], "[", "del.call"))

   seg.colors_0 = rep("black", length(called_segs_0) )
   seg.colors_0[called_segs_0] = "dodgerblue" 

   seg.colors_1 = rep("black", length(called_segs_1) )
   seg.colors_1[called_segs_1] = "dodgerblue" 

#   del.segtabs = SCNA_calls[["del.segtabs"]]
   del_foc_0 =  unlist(lapply( called.segtabs[wgd==0], "[", "del_foc" ))
   del_foc_1 =  unlist(lapply( called.segtabs[wgd>0], "[", "del_foc" ))

   plot( cn_0, del_foc_0, main="", xlim=c(0, 4), xlab="corrected CN", ylab="Focality", pch=".", col=seg.colors_0 )
   abline( v=0, lty=3, col=2 )

   plot( cn_1, del_foc_1, main="", xlim=c(0, 4), xlab="corrected CN", ylab="Focality", pch=".", col=seg.colors_1 )
   abline( v=0, lty=3, col=2 )

   plot( rcn_0, del_foc_0, main="", xlim=c(0, 4), xlab="rescaled CN", ylab="Focality", pch=".", col=seg.colors_0 )
   abline( v=0, lty=3, col=2 )

   plot( rcn_1, del_foc_1, main="", xlim=c(0, 4), xlab="rescaled CN", ylab="Focality", pch=".", col=seg.colors_1 )
   abline( v=0, lty=3, col=2 )

   dev.off()







## gene-level plots

   gene_SCNA_dat = SCNA_calls[["SCNA_event_dat"]]
   amp.gene.data = gene_SCNA_dat[["amp.gene.data"]]
   del.gene.data = gene_SCNA_dat[["del.gene.data"]]



   pdf("drivergene_cn_vs_foc.pdf", 12, 12 )
   par(mfrow=c(2,2))
   par(  bty="n", las=1  )

## cut down to driver genes only
#   print("Driver genes missing:")
#   print( setdiff(drivers, genes) )
#   drivers = intersect(drivers, genes)

## subset to specified samples - make 2 plots
   plot_genes_SCNA_dat = function( gene.data, samples, main )
   {
      gene.data = gene.data[, samples, ]

      amp.call = as.vector(  (gene.data[,,"amp.call"]) )
      amp.cn = as.vector( (gene.data[,,"corrected_total_cn"]) )

      
      amp.relative.cn = as.vector( gene.data[,,"corrected_total_cn"] / matrix(ploidy[samples], nrow=dim(gene.data)[1], ncol=dim(gene.data)[2], byrow=TRUE)  )


      amp.foc =  as.vector( (gene.data[,,"amp_foc"]) )
      amp.len.kb =  as.vector( (gene.data[,,"length"]) ) / 1e3

      gene_color = rep( "black", length(amp.call) )
      gene_color[as.logical(amp.call)] = "red"

#   plot( log(amp.cn, 2), amp.foc, main="", xlim=c(0, max(log(amp.cn,2))),  xlab="Log2 corrected CN", ylab="Focality", pch=".", col=col)
      plot( 0, type="n", main=main, xlim=c(0, max(log(amp.cn,2))), ylim=c(0.975,1), xlab="Log2 corrected CN", ylab="Focality")
      text( x=log(amp.cn, 2), y=amp.foc, labels=gene_list, font=3, cex=0.5, col=gene_color)

      plot( 0, type="n", main=main, xlim=c(0, max(log(amp.relative.cn,2))), ylim=c(0.975,1), xlab="Log2 (corrected CN / ploidy)", ylab="Focality")
      text( x=log(amp.relative.cn, 2), y=amp.foc, labels=gene_list, font=3, cex=0.5, col=gene_color)

#      plot( 0, type="n", main=main, xlim=c(0, max(log(amp.cn,2))), ylim=c(0,5000), xlab="Log2 corrected CN", ylab="Length (kb)")
#      text( x=log(amp.cn, 2), y=amp.len.kb, labels=gene_list, font=3, cex=0.5, col=gene_color)
   }  


# cut down to freq amps
   gene_amp_freq = rowSums( amp.gene.data[,,"amp.call"] )
   amp_genes = names(gene_amp_freq)[gene_amp_freq >= 5]
   gene_list = intersect( amp_genes, drivers) 

   driver.amp.gene.data = amp.gene.data[gene_list,,]

   wgd_1_snames = dimnames(amp.gene.data)[[2]][wgd==1]
   wgd_0_snames = dimnames(amp.gene.data)[[2]][wgd==0]

   plot_genes_SCNA_dat( driver.amp.gene.data, wgd_0_snames, "Amps in non-doubled genomes"  )
   plot_genes_SCNA_dat( driver.amp.gene.data,  wgd_1_snames, "Amps in doubled (1) genomes"  )


## Output 1 table for each del gene:   del.dat X sample
   out.dir = ("CN/genetabs")
   dir.create(out.dir, recursive = TRUE)
   for( i in 1:length(gene_list) )
   {
      fn = file.path( out.dir, paste( gene_list[i], ".amp.table.txt", sep="") )
      gd = driver.amp.gene.data[gene_list[i],,]
      ix= order(gd[,"amp.call"], gd[,"corrected_total_cn"], gd[,"amp_foc"])
      write.table( gd[ix,], file=fn, sep="\t", quote=FALSE )
   }


#   dev.off()


## Deletions
## cut down to driver genes only
#   print("Driver genes missing:")
#   print( setdiff(drivers, genes) )
#   drivers = intersect(drivers, genes)

# cut down to freq dels
   gene_del_freq = rowSums( del.gene.data[,,"del.call"] )
   del_genes = names(gene_del_freq)[gene_del_freq >= 5]
   gene_list = intersect( del_genes, drivers) 
#   gene_list = grep("TTT", del_genes, invert=TRUE, value=TRUE)
   



   driver.del.gene.data = del.gene.data[gene_list,,]
   del.call = as.vector(  (driver.del.gene.data[,,"del.call"]) )
   del.cn = as.vector( (driver.del.gene.data[,,"corrected_total_cn"]) )
   del.rcn = as.vector( (driver.del.gene.data[,,"rescaled_total_cn"]) )
   del.foc =  as.vector( (driver.del.gene.data[,,"del_foc"]) )
   del.len.kb =  as.vector( (driver.del.gene.data[,,"length"]) ) / 1e3

   gene_color = rep( "black", length(del.call) )
   gene_color[as.logical(del.call)] = "dodgerblue"

   plot( 0, type="n", main="", xlim=c(-0.25, 1.25), ylim=c(0.975, 1), xlab="Corrected CN", ylab="Focality")
   text( x=del.cn, y=del.foc, labels=gene_list, font=3, cex=0.5, col=gene_color)

   plot( 0, type="n", main="", xlim=c(min(del.rcn), 1.25), ylim=c(0.975, 1), xlab="Rescaled CN", ylab="Focality")
   text( x=del.rcn, y=del.foc, labels=gene_list, font=3, cex=0.5, col=gene_color)


   plot( 0, type="n", main="", xlim=c(-0.25, 1.25), ylim=c(min(del.len.kb), 500), xlab="Corrected CN", ylab="Segment length (kp)")
   text( x=del.cn, y=del.len.kb, labels=gene_list, font=3, cex=0.5, col=gene_color)

   plot( 0, type="n", main="", xlim=c(min(del.rcn), 1.25), ylim=c(min(del.len.kb), 500), xlab="Rescaled CN", ylab="Segment length (kp)")
   text( x=del.rcn, y=del.len.kb, labels=gene_list, font=3, cex=0.5, col=gene_color)
  
   
  # plot del len vs. focality
   plot( 0, type="n", main="", xlim=c(0.975, 1), ylim=c(min(del.len.kb), 500), xlab="Focality", ylab="Segment length (kp)")
   text( x=del.foc, y=del.len.kb, labels=gene_list, font=3, cex=0.5, col=gene_color)
 


   dev.off()


## Output 1 table for each del gene:   del.dat X sample
   out.dir = ("CN/genetabs")
   dir.create(out.dir)
   for( i in 1:length(gene_list) )
   {
      fn = file.path( out.dir, paste( gene_list[i], ".del.table.txt", sep="") )
      gd = driver.del.gene.data[gene_list[i],,]
      ix= order(gd[,"del.call"], gd[,"corrected_total_cn"], gd[,"del_foc"])
      write.table( gd[ix,], file=fn, sep="\t", quote=FALSE )
   }
}



# also calls ABS segs: temp solution until this is done in ABS at extraction
# produce matrix of sample X region  true/false values + metadata
# Takes direction into account - e.g. only call amps in amp regs

# NOT USED currently - qualify events at the gene level instead
genotype_amp_and_del_SCNAs_in_called_ABS_files = function( ABS_BASE_DIR, amp_regs, del_regs )
{
   files = grep(  ".ABSOLUTE.SLC.called.RData", dir(ABS_BASE_DIR, full.names=FALSE), value=TRUE)
   snames = gsub( ".ABSOLUTE.SLC.called.RData", "", files )

   amp_ev_mat = matrix( NA, nrow=length(amp_regs), ncol=length(snames) )
   del_ev_mat = matrix( NA, nrow=length(del_regs), ncol=length(snames) )

   rownames(amp_ev_mat) = names(amp_regs)
   rownames(del_ev_mat) = names(del_regs)
   colnames(amp_ev_mat) = colnames(del_ev_mat) = snames

### Debugging:
#   i = 24
   i = grep( "PB0274-TM", files )
   load( file.path(ABS_BASE_DIR, files[i]) ) 
   ABS.dat = seg.obj
   called.segtab = call_genome_wide_ABSOLUTE_SCNAs( ABS.dat )
   del.reg.segtab = call_regions_ABSOLUTE_SCNAs( called.segtab, del_regs ) [["del.reg.segtab"]]
###

   res = foreach( i = 1:length(files)) %dopar%
   {
      load( file.path(ABS_BASE_DIR, files[i]) ) 
      ABS.dat = seg.obj

      called.segtab = call_genome_wide_ABSOLUTE_SCNAs( ABS.dat )
      del.reg.segtab = call_regions_ABSOLUTE_SCNAs( called.segtab, del_regs ) [["del.reg.segtab"]]
      amp.reg.segtab = call_regions_ABSOLUTE_SCNAs( called.segtab, amp_regs ) [["amp.reg.segtab"]]

      cat(".")
      return( list("called.segtab"=called.segtab, "amp.segtab"=amp.reg.segtab, "del.segtab"=del.reg.segtab) )
   }

   called.segtabs = lapply( res, "[[", "called.segtab" )
   amp.segtabs = lapply( res, "[[", "amp.segtab" )
   del.segtabs = lapply( res, "[[", "del.segtab" )

   for( i in 1:length(files) )
   {
      amp_ev_mat[,i] = amp.segtabs[[i]][,"amp.call"]
      del_ev_mat[,i] = del.segtabs[[i]][,"del.call"]
   }
   cat("done\n")

   names(called.segtabs) = names(amp.segtabs) = names(del.segtabs) = snames
   SCNA_event_dat = list("amp_ev_mat"=amp_ev_mat, "del_ev_mat"=del_ev_mat, "amp_regs"=amp_regs, "del_regs"=del_regs )  
   
   return( list( "called.segtabs"=called.segtabs, "amp.segtabs"=amp.segtabs, "del.segtabs"=del.segtabs, "SCNA_event_dat"=SCNA_event_dat) )
}


