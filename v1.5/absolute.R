## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

RunAbsolute = function(seg.dat.fn, primary.disease, platform, sample.name, results.dir, copy_num_type, genome_build, gender=NA, min.ploidy=1, max.ploidy=8, max.as.seg.count=1500, max.non.clonal=0.8, max.neg.genome=0.005, maf.fn = NULL, indel.maf.fn = NULL, min.mut.af = NULL, output.fn.base=NULL, min_probes=10, max_sd=100, sigma.h=0.01, SSNV_skew=1, filter_segs=TRUE, force.alpha=NA, force.tau=NA, verbose = FALSE) 
{  
  # new addition
  load("/xchip/tcga/Tools/absolute/releases/v1.5/data/ChrArmsDat.RData")

  compute_cached = function( results.dir, file.base, file.ext, fn, verbose, ... )
  {
 ## Caching for mode.tab
     result_FN = file.path(results.dir, paste(file.base, "_", file.ext, ".Rds", sep = ""))
     if( file.exists(result_FN))
     {
        if( verbose ) { print( paste( "loading cached ", file.ext, " result", sep="")) }
        obj = try(readRDS( result_FN))

        if( class(obj) == "try-error" ) { print( "load failed!!"); cached=FALSE }
        else { cached = TRUE }
     }
     else{ cached = FALSE }
     if( !cached )
     {
        if( verbose ) { print( paste( "Computing ", file.ext, " result", sep="")) }
        obj = fn(..., verbose=verbose)
        saveRDS( obj, file=result_FN )
     }
     return( obj )
  }


  genome_build = match.arg(genome_build, c("mm9", "hg18", "hg19") )
  
 ## TODO:  1) build a new table for hg19
 ##        2) move data from current genome.R into hg RData files
 ##	   3) genome_HSCR_seg_plot.R is currently fixed to hg18 data (in genome.R)	
  if( genome_build %in% c("hg18") ) { load("/xchip/tcga/Tools/absolute/releases/v1.5/data/hg18_ChrArmsDat.RData") }
  if( genome_build == "mm9" ) { load("/xchip/tcga/Tools/absolute/releases/v1.5/data/mm9_ChrArmsDat") }

  platform = match.arg(platform, c("SNP_6.0", "Illumina_WES"))
  if (platform == "SNP_6.0") {
    filter_segs = TRUE
  } else if (platform == "Illumina_WES") {
    filter_segs = TRUE
  } else {
    stop("Unsupported platform: ", platform)
  }
  
  if (copy_num_type == "total") {
    set_total_funcs()
  } else if (copy_num_type == "allelic") {
    set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }

## Note - soon we will switch to ASCII HAPSEG output for 6.0, then the HAPSEG functions below will be deprecated and replaced by the Allelic versions   The code block below can then be replaced by platform_funcs.R
  if( platform == "SNP_6.0" ) 
  {
## extract segtab form .RData binary
    segtab = extract_HAPSEG_segtab(seg.dat.fn, verbose=verbose) 
    MakeSegObj <<- AllelicMakeSegObj

#    MakeSegObj <<- HAPSEGMakeSegObj
  }
  else
  {
    if (!file.exists(seg.dat.fn)) {
      stop("seg.dat.fn does not exist")
    }
    segtab = read.delim( seg.dat.fn, row.names=NULL, stringsAsFactors=FALSE, check.names=FALSE)

    if (copy_num_type == "total") {
      MakeSegObj <<- total_make_seg_obj
    } else if (copy_num_type == "allelic") {
      MakeSegObj <<- AllelicMakeSegObj
    }
  }

  ##  set up SCNA and SSNV model parameters
  SCNA.argv = list( copy_num_type, min.ploidy, max.ploidy, sigma.h )
  names(SCNA.argv) = c( "copy_num_type", "min.ploidy", "max.ploidy", "sigma.h" )
  SCNA_model = SCNA_model_setup( SCNA.argv, verbose )
  SCNA_model[["kQ"]] <- max.ploidy + 1
  
  tmp.dir = file.path(results.dir, "tmp")
  dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results.dir, recursive=TRUE, showWarnings=FALSE)
  file.base = paste(output.fn.base, ".ABSOLUTE", sep = "")

  
  seg.dat = MakeSegObj(segtab, gender, min_probes=min_probes, max_sd=max_sd, 
                        filter_segs=filter_segs, verbose=verbose)

  seg.dat[["primary.disease"]] = primary.disease
  seg.dat[["group"]] = DetermineGroup(primary.disease)
  seg.dat[["platform"]] = platform
  seg.dat[["copy_num_type"]] = copy_num_type
  seg.dat[["sample.name"]] = as.character(sample.name)
  if (is.null(seg.dat$array.name)) {
    seg.dat$array.name = seg.dat$sample.name
  }
  seg.dat[["maf.fn"]] = maf.fn
  seg.dat[["indel.maf.fn"]] = indel.maf.fn
  
  ## either allelic or total CR, to be modeled.
  seg.dat[["obs.scna"]] = ExtractSampleObs(seg.dat)
  SCNA_model[["N_probes"]] = seg.dat[["obs.scna"]][["N_probes"]]

  ## check for QC failure modes
  if (verbose) {
    print(paste("Expected copy-ratio = ", round( seg.dat[["obs.scna"]][["e.cr"]], 5), sep=""))
  }
  
  mode.res = list(mode.flag = NA)
  
  if ( length(seg.dat[["obs.scna"]][["W"]] ) > max.as.seg.count) {
    mode.res[["mode.flag"]] = "OVERSEG"
  }
  
  if ((seg.dat[["obs.scna"]][["e.cr"]] < 0.5) || (seg.dat[["obs.scna"]][["e.cr"]] > 1.5)) {
    mode.res[["mode.flag"]] = "E_CR_SCALE"
  }
  
  if (is.na(mode.res[["mode.flag"]])) {
    ## check for MAF describing somatic mutations
    maf = NULL
    if ((!is.null(maf.fn)) && (file.exists(maf.fn))) {
      maf = read.delim(maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                        check.names = FALSE, na.strings = c("NA", "---"),
                        blank.lines.skip=TRUE, comment.char="#")
    } else {
       stop(paste("MAF file: ", maf.fn, " not found.", sep = ""))
    }
    
    indel.maf = NULL
    if (!is.na(indel.maf.fn) && file.exists(indel.maf.fn)) {
      indel.maf = read.delim(indel.maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                        check.names = FALSE, na.strings = c("NA", "---"),
                        blank.lines.skip=TRUE, comment.char="#")
    } else {
       print(paste("Indel MAF file: ", indel.maf.fn, " not found.", sep = ""))
    }

## find initial purity/ploidy solutions for debugger

#    data(ChrArmsDat, package = "ABSOLUTE")
    if ((!is.null(maf)) && (nrow(maf) > 0)) 
    {
      SSNV_model = init_SSNV_model( SCNA_model[["kQ"]], SSNV_skew, nrow(maf) )
      mut.cn.dat = CreateMutCnDat(maf, indel.maf, seg.dat, min.mut.af, verbose=verbose)
    } else {
      mut.cn.dat = NA
      SSNV_model = NA
    }

  ## Caching for mode.tab
    mode.tab = compute_cached( results.dir, file.base, "mode.tab", ProvisionalModeSweep, verbose, 
               seg.dat, SCNA_model, mut.cn.dat, SSNV_model, force.alpha, force.tau, chr.arms.dat )

## Caching for mode.res
    mode.res = compute_cached( results.dir, file.base, "mode.res", fit_modes_SCNA_models, verbose, 
                               seg.dat, mode.tab, SCNA_model, mut.cn.dat )
  }

  if (is.na(mode.res[["mode.flag"]])) 
  {
    bad.ix = GenomeHetFilter(seg.dat[["obs.scna"]], mode.res, max.non.clonal,
                              max.neg.genome, SCNA_model[["kQ"]], verbose=verbose)
    if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      mode.res = list(mode.flag="ALPHA_TAU_DOM")
    } else {
      mode.res = ReorderModeRes(mode.res, !bad.ix)
    }
  }

  
  if (is.na(mode.res[["mode.flag"]]))
  {
    ## 1 - apply karyotype model
    ## Kar model only defined for human cancers
    if( genome_build %in% c("hg18", "hg19") ) 
    {
       #data(ChrArmPriorDb, package="ABSOLUTE")
       load("/xchip/tcga/Tools/absolute/releases/v1.5/data/ChrArmPriorDb.RData")
       model.id = ifelse(seg.dat[["group"]] %in% names(train.obj),
                         seg.dat[["group"]], "Primary")
       mode.res = ApplyKaryotypeModel(mode.res, model.id, train.obj, verbose=verbose)
    }
    else{ seg.dat[["group"]]="" }
 
    ## 2 - apply mutation model
    if ((!is.null(maf)) && (nrow(maf) > 0)) 
    {
      seg.dat[["mut.cn.dat"]] = mut.cn.dat

## Caching for updated mode.res with SSNV results
      mode.res = compute_cached( results.dir, file.base, "SSNV.mode.res", ApplySSNVModel, verbose, 
                                 mode.res, mut.cn.dat, SSNV_model )
    }

    mode.res = WeighSampleModes(mode.res)
    mode.res[["call.status"]] = GetCallStatus(mode.res, seg.dat[["obs.scna"]][["W"]])
  }
  
  seg.dat[["mode.res"]] = mode.res
  rm(mode.res); gc()   ## try to save some mem

#  if (is.null(output.fn.base)) {
#    output.fn.base = ifelse(is.null(seg.dat$array.name), sample.name, seg.dat$array.name)
#  }
    
  save(seg.dat, file = file.path(results.dir, paste(file.base, "RData", sep = ".")))

## plot result
  if (is.na(seg.dat[["mode.res"]][["mode.flag"]])) {
    sample.pdf.fn = file.path(results.dir,
                               paste(file.base, "plot.pdf", sep = "_"))
    if( verbose ) { print("Making result plot") }

    AbsoluteResultPlot(sample.pdf.fn, seg.dat)
  } else {
    if (verbose) {
      print("Mode flag is NA, not generating plots. Sample has failed ABSOLUTE")
    }    
  }
  
  return(TRUE)
}


