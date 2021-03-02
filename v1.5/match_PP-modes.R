## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

#edited for FH standalone



run_PP_calls_liftover_from_num = function( solution_num, analyst.id, modes.fn, 
                                  out.dir.base, obj.name, verbose=FALSE) {
  ## provides segobj.list
  if( verbose ) { cat("Loading ABS multi-sample object...") }
  load(modes.fn)  
  if( verbose ) { cat(" Done.\n") }

  review.dir = paste( out.dir.base, "/reviewed/", sep="")
  dir.create(review.dir, recursive=TRUE)
  call_override=list(solution_num)
  cat(paste("Solution num:",solution_num))
 
  called.segobj.list = override_absolute_calls(segobj.list, call_override)
  
  ## PP tab
  out.fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".", analyst.id, 
                                                     ".ABSOLUTE.table.txt", sep=""))
  PrintPpCallTable(called.segobj.list, out.fn)

  return(called.segobj.list)

 }


run_PP_calls_liftover = function( solution_num, analyst.id, modes.fn, 
                                  out.dir.base, obj.name, verbose=FALSE)
{
  ## provides segobj.list
  if( verbose ) { cat("Loading ABS multi-sample object...") }
  load(modes.fn)  
  if( verbose ) { cat(" Done.\n") }

  #dat = read.delim(reviewed.pp.calls.fn, row.names=NULL, stringsAsFactors=FALSE, header=1,
  #                  check.names=FALSE)

  review.dir = paste( out.dir.base, "/reviewed/", sep="")
  dir.create(review.dir, recursive=TRUE)

if (colnames(dat)[3] == "sample") {
     print("Running PP calls override")

     call_override = dat[, 1]
  
     pp.calls = dat[, c(3:ncol(dat))]
     rownames(pp.calls) = dat[, 2]
     names(call_override) = dat[,2]

     found = intersect(names(segobj.list), rownames(pp.calls))
     missing = setdiff( rownames(pp.calls), names(segobj.list) )
     segobj.list = segobj.list[found]
     call_override = call_override[found]

     msg = paste( length(missing), " samples in pp.calls file not found in segobj.list: ", sep="")
     print(msg)
     print(missing) 

     called.segobj.list = override_absolute_calls(segobj.list, call_override)
  } else {
     if (colnames(dat)[2] != "sample") { 
       stop("Invalid reviewed.pp.calls.fn!") 
     }

     print("Running PP calls liftover...")
     pp.calls = dat[, c(2:ncol(dat))]
     rownames(pp.calls) = dat[, 1]

     found = intersect(names(segobj.list), rownames(pp.calls))

     if( length(found) == 0 )
     {
        print("No matching PP calls to liftover, exiting.")
        return()
     }

     segobj.list = segobj.list[found]
     match.dat = MatchPpModes(pp.calls, segobj.list)
     mode.ix= match.dat[,"mode.ix"]     
     nix = !is.na(mode.ix)
     wrong.ix = which(match.dat[nix, "mode.ix"] != 1)

     called.segobj.list = ReduceSegobjListToCalled(mode.ix, segobj.list)
   
     msg = paste( "Matched ", sum(!is.na(mode.ix)), " of ", length(mode.ix), " samples", sep="")
     print(msg)
     auto_pct = round( ((sum(nix)-length(wrong.ix)) / sum(nix)) * 100, 2 )
     msg = paste( length(wrong.ix), " / ", sum(nix), " matches with alternate solution top-ranked (auto call rate = ", auto_pct, "%)", sep="")
     print(msg)

     PP_liftover_plot( match.dat, paste( review.dir, obj.name, "_PP-liftover.pdf", sep=""))
    
     print("Unmatched:")
#     print( found[ is.na(mode.ix)])
#
     nix = which(is.na(mode.ix))
     if( length(nix) > 0 )
     {
        unmatched_sid = rep(NA, length(nix))
        for( i in seq_along(nix))
        { 
           unmatched_sid[i] = segobj.list[[nix[i]]][["sample.name"]]
        } 
        print( unmatched_sid )

      ## PP tab
        unmatched.out.fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".unmatched.", analyst.id, ".ABSOLUTE.table.txt", sep="") )
        PrintPpCallTable(segobj.list[nix], unmatched.out.fn)

  ## print top 3 solutions for all samples where no match found (1 plot)
        unmatched.pdf.fn = paste( review.dir, obj.name, "_unmatched_PP-modes_plot.pdf", sep="")
        pdf( unmatched.pdf.fn, 17.5, 18.5 )
        PlotModes_layout()
        plot.wix = nix

        for( i in 1:length(plot.wix)) 
        {        
           PlotModes(segobj.list[[plot.wix[i]]], n.print=3)
        }
        dev.off()
     }
  }


 ## PP tab
  out.fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".", analyst.id, 
                                                     ".ABSOLUTE.table.txt", sep=""))
  PrintPpCallTable(called.segobj.list, out.fn)

  return(called.segobj.list)
}



override_absolute_calls = function(segobj.list, call_override, verbose=FALSE) {
  mode_str = call_override
  mode.ix = rep(NA, length(mode_str))
  status.vec = rep(NA, length(mode_str))

  empty_ix = (call_override == "") | (is.na(call_override))
  call_override[empty_ix] = "1"
  called_ix = !is.na(as.integer(call_override))
  status.vec[called_ix ] = "called"   # PP_calls[ix,"call status"]

  mode.ix[called_ix] = as.integer(call_override[called_ix])
  mode.ix[!called_ix] = 1
  status.vec[!called_ix] = call_override[!called_ix] 

  segobj.list = SelectMatchedModes(segobj.list, mode.ix, status.vec, verbose=verbose)
  
  return( segobj.list )
}






MatchPpModes <- function(pp.calls, segobj_list) {
  N <- length(segobj_list)
  mode.ix <- rep(NA, N)
  wrong.ix <- c()
  names(mode.ix) <- names(segobj_list)

  match.dat = matrix(NA, nrow=N, ncol=5)
  colnames(match.dat)=c("mode.ix", "best_ploidy", "best_purity", "called_ploidy", "called_purity")
  
  for (i in seq_len(N)) 
  {
    sid = names(segobj_list)[i]
    seg.dat = segobj_list[[i]]

    if ((any(is.na(pp.calls[sid, c("purity", "ploidy")]))) ||
      (!is.na(seg.dat[["mode.res"]][["mode.flag"]]))) 
    { next }

#    mode.ix[i] <- MatchPpMode(pp.calls, sid, seg.dat)
    match.dat[i,] <- MatchPpMode(pp.calls, sid, seg.dat)

#    if(!is.na(mode.ix[i]) && mode.ix[i]!=1)
#    {
#      wrong.ix <- c(wrong.ix, i)
#    }
  }
  
#  return(list(matched=mode.ix, wrong.ix=wrong.ix))
  return(match.dat)
}

old_MatchPpMode <- function(pp.calls, sid, seg.dat, ploidy_tol=0.25, purity_tol=0.05)
{
  mode.tab <- seg.dat[["mode.res"]][["mode.tab"]]
  vals <- mode.tab[, c("alpha", "genome mass"), drop=FALSE]   
  
  call.vals <- pp.calls[sid, c("purity", "ploidy")]
  call.vals <- matrix(as.numeric(call.vals), ncol=2,
                      nrow=nrow(mode.tab), byrow=TRUE)
  
  scales = matrix(c(25,1),nrow=nrow(vals), ncol=2, byrow=TRUE )
  min.mode <- which.min(rowSums(  scales*(vals - call.vals)^2))
  
## finds closest match and accepts within epsilon.
## TODO: a better solution might be to find all modes in epsilon neighborhood and select highest LL

  best_ploidy = vals[min.mode, "genome mass"]
  best_purity = vals[min.mode, "alpha"]

  if ((abs(pp.calls[sid, "ploidy"] - best_ploidy) > ploidy_tol) |
      (abs(pp.calls[sid, "purity"] - best_purity) > purity_tol))
  {
    mode.ix=NA
  }
  else{ mode.ix <- min.mode 
  }

  return( c(mode.ix,best_ploidy,best_purity,pp.calls[sid, "ploidy"], pp.calls[sid, "purity"]  ) )
}

MatchPpMode <- function(pp.calls, sid, seg.dat, ploidy_tol=0.25, purity_tol=0.05)
{
  mode.tab <- seg.dat[["mode.res"]][["mode.tab"]]
  vals <- mode.tab[, c("alpha", "genome mass"), drop=FALSE]   
  
  eps.N.ix = (abs(pp.calls[sid, "ploidy"] - vals[,"genome mass"]) <= ploidy_tol) &
             (abs(pp.calls[sid, "purity"] - vals[,"alpha"]) <= purity_tol)

  if( sum( eps.N.ix) > 0 )
  {
    mode.ix = which(eps.N.ix)[ which.max( mode.tab[eps.N.ix, "combined_LL"] ) ]
    best_ploidy = mode.tab[mode.ix,"genome mass"]
    best_purity = mode.tab[mode.ix,"alpha"]
  }
  else
  {
    mode.ix=NA

    call.vals <- pp.calls[sid, c("purity", "ploidy")]
    call.vals <- matrix(as.numeric(call.vals), ncol=2, nrow=nrow(mode.tab), byrow=TRUE)
    scales = matrix(c(25,1),nrow=nrow(vals), ncol=2, byrow=TRUE )
    min.mode <- which.min(rowSums(  scales*(vals - call.vals)^2))
  
    best_ploidy = vals[min.mode, "genome mass"]
    best_purity = vals[min.mode, "alpha"]
  }

  return( c(mode.ix, best_ploidy, best_purity, pp.calls[sid, "ploidy"], pp.calls[sid, "purity"]  ) )
}


ReduceSegobjListToCalled <- function(mode.ix, segobj_list) {
  segobj_list <- segobj_list[!is.na(mode.ix)]
  mode.ix <- mode.ix[!is.na(mode.ix)]
  
  for (i in seq_along(segobj_list))  {
    segobj_list[[i]][["mode.res"]] <- ReorderModeRes(segobj_list[[i]][["mode.res"]],
                                                     mode.ix[i])
  }
  
  return(segobj_list)
}

SelectMatchedModes <- function(segobj_list, mode.ix, status.vec, verbose=FALSE) {
   new.segobj.list <- segobj_list[!is.na(mode.ix)]
   mode.ix <- mode.ix[!is.na(mode.ix)]

   for (i in seq_along(mode.ix)) {
      segobj <- new.segobj.list[[i]]
      segobj[["mode.res"]][["call.status"]] <- status.vec[i]

      if (!(status.vec[i] %in% c("non-clonal", "non-aneuploid",
                                "low purity", "FAILED"))) {
         mode.res <- segobj[["mode.res"]]
         ix <- mode.ix[i]
         segobj[["mode.res"]] <- ReorderModeRes(mode.res, ix)
      }

      new.segobj.list[[i]] <- segobj
      if (verbose) {
        cat(".")
      }
   }
   if (verbose) {
     cat("\n")
   }
   
   return(new.segobj.list)
}
