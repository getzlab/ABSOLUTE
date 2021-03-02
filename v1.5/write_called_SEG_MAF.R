## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

write_called_seg_maf = function(called_segobj_list, pp_calls, out_dir, verbose=FALSE) {
  dir.create(out_dir) 
  
  for (i in seq_along(called_segobj_list)) {
    called_segobj = called_segobj_list[[i]] 
    s_name = called_segobj$array.name
    
    abs_seg = GetAbsSegDat(called_segobj)
    WriteAbsSegtab(list(abs_seg), s_name, 
                   file.path(out_dir, paste(s_name, "segtab.txt", sep=".")))
    ## MAF
    if (!is.null(called_segobj$mode.res$modeled.muts)) {
      maf_out_fn = file.path(out_dir, paste(s_name, "_ABS_MAF.txt", sep=""))
            
      modeled = called_segobj$mode.res$modeled.muts[[1]]
      mut_dat = cbind(called_segobj$mut.cn.dat, modeled)      

      SSNV_ccf_dens = called_segobj[["mode.res"]][["SSNV.ccf.dens"]][1,,]
      out_mut_dat = cbind(mut_dat, SSNV_ccf_dens)
            
      write.table(file=maf_out_fn, out_mut_dat, row.names=FALSE, sep="\t", quote=FALSE)
    }
    if (verbose) {   
      cat(".")
    }
  }
}

WriteAbsSegtab <- function(seg, s_name, out_fn) {
  for (s in seq_along(seg)) {
    sample <- rep(s_name, nrow(seg[[s]]))
    s_tab <- cbind(sample, seg[[s]])
    
    ## colames only for 1st sample
    app <- s > 1   
    col <- s == 1
    
    write.table(s_tab, file=out_fn, col.names=col, append=app, row.names=FALSE,
                quote=FALSE, sep="\t")
  }   
}

