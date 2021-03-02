## functions for running ABSOLUTE with total CR

total_make_seg_obj = function(dat_fn, gender, filter_segs=FALSE, min_probes=NA, max_sd=NA, verbose=FALSE) {   
  segs_tab = read.delim(dat_fn, row.names=NULL, check.names=FALSE, stringsAsFactors=FALSE)
  seg_dat = list()
  
  print("Gender not supported yet!!")
  gender = NA

  if( !is.na(gender) && gender %in% c( "Male", "Female") )
  {
     nix = segs_tab[,"Chromosome"] %in% c("Y", "M", "chrY", "chrM")
  }
  else
  {
     if( verbose ) {
       print("No or invalid gender specified - dropping X chromosome segs")
     }
     nix = segs_tab[,"Chromosome"] %in% c("X", "Y", "M", "chrX", "chrY", "chrM")
  }
  seg_dat$gender = gender
  segs_tab = segs_tab[!nix, ]



 ## NEED this for mm9 
  segs_tab[,"Chromosome"] = as.integer(gsub( "chr", "", segs_tab[,"Chromosome"]))
  if( any(is.na(segs_tab[,"Chromosome"]))) { stop("converted to non-integer chromosome") }



  segtab = segs_tab[, c("Chromosome", "Start", "End", "Num_Probes")]
  colnames(segtab) = c("Chromosome", "Start.bp", "End.bp", "n_probes")
  
  length = segtab[, "End.bp"] - segtab[, "Start.bp"]
  ## Convert from base 2 log
  copy_num = 2^(segs_tab[, "Segment_Mean"] )   
  
  ix = copy_num > 5.0
  if (verbose) {
    print( paste( "Capping ", sum(ix), " segs at tCR = 5.0", sep=""))
  }
  copy_num[ix] = 5.0
  
  seg_sigma_num = 0.1  ## TODO - get rid of this - not used in downstream model - but crashes filtering code if missing
  seg_sigma =  seg_sigma_num / sqrt(as.numeric(segs_tab[,"Num_Probes"]))   
#  seg_sigma = rep(NA, nrow(segs_tab))  ## calculate later in SCNA_model
  
  segtab = cbind(segtab, length, copy_num, seg_sigma )

  if (filter_segs) {
    seg_dat$segtab = FilterSegs(segtab, min_probes=min_probes, max_sd=max_sd)$seg.info
  }

  W = as.numeric(seg_dat$segtab[,"length"])
  W = W / sum(W)
  seg_dat$segtab = cbind(seg_dat$segtab, W)
  colnames(seg_dat$segtab)[ncol(seg_dat$segtab)] = "W"

  seg_dat$error_model = list()

  return(seg_dat)  
}

total_extract_sample_obs = function(seg.dat) {

  seg = seg.dat[["segtab"]]
  d = seg[, "copy_num"]

  ## expected copy-number, should be 1.0
  e.cr = sum(seg[, "W"] * d )

  if( !is.na(seg.dat$gender) && seg.dat$gender == "Male" ) 
  {
     male_X = (seg[,"Chromosome"]=="X")
  }
  else
  {
     male_X = rep(FALSE, length(d) )
  }

  ## FIXME: "error.model" was originally named "HSCN_params" - double check this
  obs = list(d=d, d.tx=d, W=seg[,"W"], n_probes=seg[,"n_probes"], seg.ix=seq_along(d), error.model=seg.dat$error_model, e.cr=e.cr, data.type="TOTAL", platform=seg.dat[["platform"]], male_X=male_X )
  
  return(obs)
}

CAPSEG_get_seg_sigma = function(SCNA_model, obs)
{
  seg_sigma =  exp(SCNA_model[["Theta"]]["sigma.A"]) / (sqrt(obs[["n_probes"]] ))   
  sigma.h = SCNA_model[["sigma.h"]]

  seg_sigma = sqrt(seg_sigma^2 + sigma.h^2)
}
