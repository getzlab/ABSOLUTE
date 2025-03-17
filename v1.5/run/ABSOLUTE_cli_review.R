library(optparse)
option.list <- list(
  make_option("--file_name", dest="file_name"),
  make_option("--sample_name", dest="sample_name"),
  make_option("--results_dir", dest="results_dir", default=getwd()),
  make_option("--max_ploidy", dest="max_ploidy",type="double",default=6.0)
)

#--seg_dat_fn --maf_fn --indelmaf_fn --sample_name --results_dir --ssnv_skew

opt <- parse_args(OptionParser(option_list=option.list))
print(opt)

min.ploidy= 1.1
max.ploidy= opt[["max_ploidy"]]
sigma.h= 0.01
max.non.clonal= 0.99
min_probes= 1
max.as.seg.count= 5000
primary.disease= NA
platform= "Illumina_WES"
copy_num_type= "allelic"
genome_build= "hg19"
#CGA_DIR= "~/CGA/R/"
max.neg.genome= 0.05
maf.fn= opt[["maf_fn"]]
indel.maf.fn= opt[["indelmaf_fn"]]
min.mut.af= 0
output.fn.base= opt[["sample_name"]]
max_sd= 100
SSNV_skew= opt[["ssnv_skew"]]
filter_segs= TRUE
force.alpha= NA
force.tau= NA
verbose= TRUE
results.dir= opt[["results_dir"]]
sample.name= opt[["sample_name"]]
file_name= opt[["file_name"]]
gender= NA
#library(ABSOLUTE)


require(GenomicRanges)
require(gplots)
require(RColorBrewer)

#suppressPackageStartupMessages(require(gplots))
#suppressPackageStartupMessages(require(RColorBrewer))
CGA_DIR_ABS="/xchip/cga_home/igleshch/soft/local/absolute"

   print( paste("sourcing files in ", CGA_DIR_ABS, sep=""))
   #rr = dir( file.path( CGA_DIR_ABS, "ABSOLUTE/sandbox"), full.names=TRUE )
   rr = dir(CGA_DIR_ABS,full.names=TRUE )
   for( i in 1:length(rr) ) { 
   source(rr[i])}



   source( file.path( CGA_DIR_ABS, "Phylogic", "log_DP_post_fixed.R" ) )
   source( file.path( CGA_DIR_ABS, "Phylogic", "log_DP_post.R" ) )
   source( file.path( CGA_DIR_ABS, "Phylogic", "run_DP.R" ) )
   source( file.path( CGA_DIR_ABS, "Phylogic", "DP_utils.R" ) )
   source( file.path( CGA_DIR_ABS, "Phylogic", "summarize_DP.R" ) )
   source( file.path( CGA_DIR_ABS, "Phylogic", "CCF_DP_plots.R" ) )



#debug (RunAbsolute)

#### Step 2: collect results for summarization:
file.base = paste(file.path(results.dir,sample.name, full.names=TRUE ) , ".ABSOLUTE", sep = "")
absolute.files= file.path( results.dir, paste(file.base, "RData", sep = "."))
print(absolute.files)

CreateReviewObject(file_name , absolute.files, results.dir, "allelic", plot.modes=TRUE, num_solutions_plotted=3, verbose=TRUE) 
#####

}
 

apply_review_and_extract = function( pp.review.fn, obj.name, analyst.id )
{
  out.dir.base = file.path( "ABSOLUTE_results", obj.name )

  nm = file.path(out.dir.base, paste(obj.name, ".PP-modes", sep = ""))
  modesegs.fn = paste(nm, ".data.RData", sep = "")

  if( file.exists(pp.review.fn) )
  {
    ExtractReviewedResults( pp.review.fn, analyst.id, modesegs.fn, out.dir.base, obj.name, verbose=FALSE )
  }
}


set_allelic_funcs()
calls_FN = file_name
#debug(apply_review_and_extract)
apply_review_and_extract( pp.review.fn=calls_FN, obj.name=sample.name, analyst.id="IL" )
