## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


## Note: this code should be deprecated by standardizing the output of HAPSEG (SNP6) to an
## allelic segtab TSV format.   Once this is done, make_allelic_seg_obj.R should replace this## file.


extract_HAPSEG_segtab <- function(seg.dat.fn, verbose=FALSE) 
{   
 AllelicInvAtten <- function(r, at) {
    return(r / (1 + at - (at * r)))
  }

  ConvertNameStyle = function(names) {
    return(tolower(gsub("_", ".", names)))
  }

  RenameSegDatFields = function(seg.dat) 
  {
  ## FIXME: Total hack, used to go back and forth between the old
  ## hapseg & current absolute
  
    names(seg.dat) = ConvertNameStyle(names(seg.dat))
    names(seg.dat[["em.res"]]) = ConvertNameStyle(names(seg.dat[["em.res"]]))
 # names(seg.dat[["obs.scna"]]) = ConvertNameStyle(names(seg.dat[["obs.scna"]]))
 # names(seg.dat[["error.model"]]) = ConvertNameStyle(names(seg.dat[["error.model"]]))
  #names(seg.dat[["mode.res"]]) = ConvertNameStyle(names(seg.dat[["mode.res"]]))
  
    return(seg.dat)
  }

  ## provides seg.dat
  MAX_TRIES = 5
  tries = 1
#  while( tries < MAX_TRIES & !file.exists(seg.dat.fn)) {
#    print( paste( "seg.dat.fn: ", seg.dat.fn, " does not exist!", sep=""))
#    tries = tries + 1
#  }

#  load(seg.dat.fn)

  while( tries < MAX_TRIES ) 
  {
    res = try( load(seg.dat.fn) )
    if( class(res)!="try-error" ) { break }
#    print( paste( "seg.dat.fn: ", seg.dat.fn, " does not exist!", sep=""))
    tries = tries + 1
  }


  if (!exists("seg.dat")) {
    stop("Invalid object contained in seg.dat.fn")
  }

  ## Detect pre-1.0 HAPSEG inputs and convert
  if ("EM_res" %in% names(seg.dat)) {
    if (verbose) {
      print("Detected a pre-1.0 HAPSEG input, converting to 1.0 format ....")
    }
    seg.dat = RenameSegDatFields(seg.dat)
  }


## Allelic CAPSEG colnames....
#Chromosome	Start.bp	End.bp	n_probes	length	f	tau	sigma.tau	mu.minor	sigma.minor	 mu.major	sigma.major	Interest
  
#colnames(seg.dat$allele.segs)
# "Chromosome" "Start.bp"   "End.bp"  "n_probes"   "length"  "A1.Seg.CN"  "A2.Seg.CN"  "tCN.Seg.sd" "AS.Seg.sd" 
  new.seg.dat <- seg.dat[["allele.segs"]]

# rescale SD
  new.seg.dat[,"AS.Seg.sd"] = sqrt( new.seg.dat[,"AS.Seg.sd"]^2 + new.seg.dat[,"tCN.Seg.sd"]^2 ) / 2

# remove attenuation and compute tCR
  AT = seg.dat[["em.res"]][["theta"]][["at"]]  
  new.seg.dat[,"A1.Seg.CN"] = AllelicInvAtten(new.seg.dat[,"A1.Seg.CN"], AT)
  new.seg.dat[,"A2.Seg.CN"] = AllelicInvAtten(new.seg.dat[,"A2.Seg.CN"], AT)
## compute total CR.  This is exactly equivalent to: AllelicInvAtten( e.mu[,3], AT )
  new.seg.dat = cbind( new.seg.dat, "total.Seg.CN" = rowSums( new.seg.dat[,c("A2.Seg.CN", "A1.Seg.CN")] ) )

## normalize column names to be equivelent to AllelicCapseg...
  segtab = data.frame( new.seg.dat[, c("Chromosome", "Start.bp", "End.bp", "n_probes", "length")], "f"=NA, 
                  "tau" = new.seg.dat[,"total.Seg.CN"],
                  "sigma.tau"=new.seg.dat[,"tCN.Seg.sd"], 
                  "mu.minor" =new.seg.dat[,"A1.Seg.CN"],
                  "sigma.minor" = new.seg.dat[,"AS.Seg.sd"],
                  "mu.major" =new.seg.dat[,"A2.Seg.CN"],
                  "sigma.major" = new.seg.dat[,"AS.Seg.sd"] )

# rows with invalid allelic CN results from HAPSEG
  HSCN_cols = c("sigma.minor", "sigma.major", "mu.minor", "mu.major")
  n.ix = apply( is.na(segtab[,HSCN_cols]), 1, any )
  segtab[ n.ix, HSCN_cols] = NA


  return( segtab )
}

