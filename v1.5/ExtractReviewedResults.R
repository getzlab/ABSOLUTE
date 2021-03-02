## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

ExtractReviewedResults = function( called.segobj.list, analyst.id, out.dir.base, obj.name, verbose=FALSE) 
{
 ## agg MAF
  cat("Outputting aggregate MAF...")
  MAF_list_fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".MAF_list.Rds", sep=""))
  if( !file.exists(MAF_list_fn) )
  {
     MAF_list = get_MAF_list_from_called_seglist_obj( called.segobj.list )
     saveRDS( MAF_list, file=MAF_list_fn )
  } else{ MAF_list = readRDS(MAF_list_fn) }

  fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".aggregate_MAF.Rds", sep=""))
  if( !file.exists(fn) )
  {
     AGG_MAF = aggregate_sample_MAF_list( MAF_list )
     saveRDS( AGG_MAF, file=fn)
  } else{ AGG_MAF = readRDS(fn) }
  cat("done\n")
 ##

 ## significance of SSNVs
#  cat("Calculating significance of SSNVs...")
#  out.base = paste( file.path(out.dir.base, "reviewed", obj.name), "_", sep="")
#  ABS_class_mutsig( AGG_MAF, called.segobj.list, out.base )
#  cat("done\n")
 ##


 ## SEG_MAFs
  cat("Extracting SEG_MAF files...")
  seg.maf.dir = file.path(out.dir.base, "reviewed", "SEG_MAF")
  dir.create(seg.maf.dir, recursive=TRUE)
  write_called_seg_maf(called.segobj.list, pp.calls, seg.maf.dir)
  cat("done\n")
 ##

 ## Called summary plot
  pdf.fn = file.path(out.dir.base, "reviewed", 
                     paste(obj.name, ".called.ABSOLUTE.plots.pdf", sep=""))
  
#  PlotModes(called.segobj.list, pdf.fn, n.print=1)

## This is useless because it does not include sample names in the plot!
if( FALSE )
{
  cat("Plotting called mode for matched samples")
  pdf( pdf.fn, 17.5, 18.5 )
  PlotModes_layout()
  for( i in 1:length(called.segobj.list)) 
  {        
     PlotModes(called.segobj.list[[i]], n.print=1)
     cat(".")
  }
  dev.off()
  cat("done\n") 
}

  ## Called indv. RData files
   cat("Extracting RData called mode files for matched samples")
   indv.called.dir = file.path(out.dir.base, "reviewed", "samples")
   dir.create(indv.called.dir, recursive=TRUE)

   file.base = file.path(paste(names(called.segobj.list), ".ABSOLUTE.", analyst.id, 
                               ".called", sep = ""))
   called.files= file.path(indv.called.dir, paste(file.base, "RData", sep = "."))
   for (i in seq_along(called.files)) {
      seg.obj = called.segobj.list[[i]]
      save(seg.obj, file=called.files[i])
      cat(".")
   }
   cat("done\n")
 
 ## print detailed SSNV plots for called solutions (1 plot per sample)
## These are intended to help debug SSNV on subclonal SCNA analysis
#   called_detailed_SSNV_plots( called.segobj.list, out.dir.base )
}

apply_review_and_extract = function( pp.review.fn, obj.name, analyst.id, verbose=TRUE )
{
   calls_FN = pp.review.fn
   out.dir.base = file.path( "ABSOLUTE_results", obj.name )

   called.out.fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".", analyst.id, 
                                                     ".called.segobj.list.Rds", sep=""))

   nm = file.path(out.dir.base, paste(obj.name, ".PP-modes", sep = ""))
   modesegs.fn = paste(nm, ".data.RData", sep = "")

   if( file.exists(calls_FN) )
   {
      called.segobj.list = run_PP_calls_liftover(calls_FN, analyst.id, modesegs.fn, out.dir.base, obj.name, verbose=verbose )
   }
   else{ stop("calls_FN does not exist!") }

   if( length(called.segobj.list) > 0 )
   {
      ExtractReviewedResults( called.segobj.list, analyst.id , out.dir.base, obj.name, verbose=TRUE )
   }
   else{ stop("called.segobj.list has length 0!") }
}


called_detailed_SSNV_plots = function( called.segobj.list, out.dir.base )
{
  plot_dir = file.path(out.dir.base, "reviewed", "SSNV_detail")
  dir.create( plot_dir, recursive=TRUE)

  for( i in 1:length(called.segobj.list) )
  {
     SID = names(called.segobj.list)[i]
     pdf.fn = file.path( plot_dir, paste(SID, ".SSNV.detail.plot.pdf", sep=""))

     seg.dat = called.segobj.list[[i]]
     mut.cn.dat <- seg.dat[["mut.cn.dat"]]
     modeled <- seg.dat[["mode.res"]][["modeled.muts"]][[1]]
     modeled.mut.dat <- cbind(mut.cn.dat, modeled)
#     SSNV_skew = modeled.mut.dat[1,"SSNV_skew"]
     SSNV_model = called.segobj.list[[i]][["mode.res"]][["mode_SSNV_models"]][[1]]

     pdf( pdf.fn, 10, 10 )

     detailed_SSNV_on_subclonal_HSCN_plot( SSNV_model, seg.dat, modeled.mut.dat, mode.ix=1, verbose=TRUE )
     
     dev.off()
  }
}


PP_liftover_plot = function( match.dat, PDF_OUT_FN )
{
   pdf( PDF_OUT_FN, 8, 4 )
   par(mfrow=c(1,2))
   par(bty='n')
   par(las=1)

   N = nrow(match.dat)
   PCH= rep(NA, N)
   col= rep(NA, N)

   nix = is.na(match.dat[,"mode.ix"])
   PCH=16
   col[nix] ="red"
   col[!nix]="black"

   CEX=0.8
   

   plot( match.dat[,"called_purity"], match.dat[,"best_purity"], xlim=c(0,1), ylim=c(0,1), main="", xlab="Called purity", ylab="Best match", pch=PCH, col=col )
   abline( a=0, b=1, col="grey", cex=CEX)


   MAXP=max(match.dat[,c("called_ploidy","best_ploidy")], na.rm=TRUE)
   plot( match.dat[,"called_ploidy"], match.dat[,"best_ploidy"], xlim=c(1,MAXP), ylim=c(1,MAXP), main="", xlab="Called ploidy", ylab="Best match",  pch=PCH, col=col, cex=CEX )
   abline( a=0, b=1, col="grey")

   dev.off()
}

