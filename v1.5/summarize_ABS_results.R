## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CreateReviewObject = function(obj.name, absolute.files, indv.results.dir, copy_num_type, plot.modes=TRUE, num_solutions_plotted=3, verbose=FALSE)
{
  
  if (copy_num_type == "total") {
    set_total_funcs()
  } else if (copy_num_type == "allelic") {
    set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  
  dir.create(indv.results.dir, recursive=TRUE)
  nm = file.path(indv.results.dir, paste(obj.name, ".PP-modes", sep = ""))
  modesegs.fn = paste(nm, ".data.RData", sep = "")
  pdf.fn = paste(nm, ".plots.pdf", sep = "")
  failed.pdf.fn = paste(nm, "FAILED_plots.pdf", sep = ".")
  failed.tab.fn = paste(nm, "FAILED_tab.txt", sep = ".")
  call.tab.fn = file.path(indv.results.dir,
                           paste(obj.name, "PP-calls_tab.txt", sep = "."))
  

  ## read in processed SEG / MODES results and assemble
  segobj.list = vector(mode = "list", length = length(absolute.files))
  failed.list = vector(mode="list", length=length(absolute.files))
  so.ix = 1
  fa.ix = 1
  
  for (i in seq_along(absolute.files)) {
    ## read in absolute rda
    seg.out.fn = absolute.files[i]
    if (!file.exists(seg.out.fn)) {
      if (verbose) {
        cat("\n")
        print(paste("sample #", i, " result not found", sep = ""))
      }
      next
    }

    ## provides seg.dat
    load(absolute.files[i])
    if (is.null(seg.dat$array.name)) {
      seg.dat$array.name = seg.dat$sample.name
    }
    SID = seg.dat$array.name
   
    if (is.na(seg.dat[["mode.res"]][["mode.flag"]])) {
      segobj.list[[so.ix]] = seg.dat
      names(segobj.list)[so.ix] = SID
      so.ix = so.ix + 1
      
      if (verbose) {
        cat(".")
      }
    } else {
      failed.list[[fa.ix]] = seg.dat
      names(failed.list)[fa.ix] = SID
      failed.list[[fa.ix]][["sample.name"]] = seg.dat[["sample.name"]]
      fa.ix = fa.ix + 1
      if (verbose) {
        cat("-")
      }
    }
  }
  
  segobj.list = segobj.list[c(1:(so.ix - 1))]
  failed.list = failed.list[c(1:(fa.ix - 1))]
  
  ## sort samples by the entropy of the best solution
  mode.ent = rep(NA, length(segobj.list))
  names(mode.ent) = names(segobj.list)
  for (i in seq_along(segobj.list)) {
    mtab = segobj.list[[i]][["mode.res"]][["mode.tab"]]

if( any( is.na(mtab[, "combined_LL"]) ) ) 
{
   cat("Bad sample found: ")
   cat(segobj.list[[i]][["sample.name"]] )
   cmd = paste( "echo ", segobj.list[[i]][["sample.name"]], " >> bad_samples", sep="")
   system(cmd)
   cat("\n")
   next
}

    ix = which.max(mtab[, "combined_LL"])
    mode.ent[i] = mtab[ix, "entropy"]
  }
  samples = names(sort(mode.ent))
  segobj.list = segobj.list[samples]
  save(segobj.list, file = modesegs.fn)


  PrintPpCallTable(segobj.list, call.tab.fn)
  
  if (plot.modes) 
#  if( FALSE )
  {
    pdf(pdf.fn, 17.5, 18.5 )
    PlotModes_layout()
    for( i in 1:length(segobj.list) )
    {
       PlotModes( segobj.list[[i]], n.print=3 )
    }
    dev.off()
  }

  if (!is.null(failed.list[[1]])) 
  {
    try( PlotFailedSamples(failed.list, failed.pdf.fn) )
    PrintFailedTable(failed.list, failed.tab.fn)
  }

  rm(segobj.list);  gc() 
  
  return(TRUE)
}
