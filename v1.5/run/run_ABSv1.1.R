library(optparse)
option.list <- list(
  make_option("--sif", dest="sif"),
  make_option("--name", dest="name"), 
  make_option("--norun", dest="dry_run",action="store_true")) 

#--seg_dat_fn --maf_fn --indelmaf_fn --sample_name --results_dir --ssnv_skew

opt <- parse_args(OptionParser(option_list=option.list))
print(opt)


source("/xchip/cga_home/igleshch/CancerGenomeAnalysis/trunk/R/broad_utils/batch_exec_ABSOLUTE.R")
source("/xchip/cga_home/igleshch/CancerGenomeAnalysis/trunk/R/broad_utils/CN_MAF_file_lookup.R")

options(error=recover)

#library(ABSOLUTE)
CGA_DIR = "/xchip/cga_home/igleshch/CancerGenomeAnalysis/trunk/R" #"/xchip/cga_home/igleshch/soft/local/absolute"
if( "CGA_DIR" %in% ls() )
{
   rr = dir( file.path( CGA_DIR, "ABSOLUTE/sandbox"), full.names=TRUE )
   for( i in 1:length(rr) ) { source(rr[i]) }
}
source( file.path( CGA_DIR, "Phylogic", "run_DP.R" ) )
source( file.path( CGA_DIR, "Phylogic", "CCF_DP_plots.R") )


### setup arguments to RunAbsolute
   ABSOLUTE_argv = list()

   ABSOLUTE_argv$min.ploidy = 1.1
   ABSOLUTE_argv$max.ploidy = 6.0
   ABSOLUTE_argv$sigma.h = 0.01
   ABSOLUTE_argv$max.non.clonal = 0.99
   ABSOLUTE_argv$min_probes = 1
   ABSOLUTE_argv$max.as.seg.count=5000

   ABSOLUTE_argv$primary.disease = NA
   ABSOLUTE_argv$platform = "Illumina_WES"
   ABSOLUTE_argv$copy_num_type = "allelic"
   ABSOLUTE_argv$genome_build = "hg19"

# not passed to RunAbsolute
   ABSOLUTE_argv$CGA_DIR = CGA_DIR
###

#fiss annot_get An_Engelman_ALKResistance pair pset=PR_Engelman_ALKResistance_PAIR_All_Pairs alleliccapseg_tsv maf_file_capture_indel_forcecalled maf_union_forcecalls AllelicCapseg_skew > FH_SIF_04.07.txt

   SIF_FN = opt[["sif"]] 
   SIF = read.delim(SIF_FN, row.names=1, check.names=FALSE, stringsAsFactors=FALSE )

   calls_FN = paste(opt[["name"]],".table.txt")
   var_bsub_argv = firehose_CAPSEG_SIF(  SIF, calls_FN, FORCE_CALL=FALSE, EXCLUDE_CALLED=FALSE, EXCLUDE_PASSED=FALSE )

   obj.name = opt[["name"]]


## Step 1: run ABSOLUTE on each sample - create a review table / plot when finished
if( TRUE )
{
   batch_exec_ABSOLUTE( ABSOLUTE_argv, obj.name, var_bsub_argv, queue="hour", dry_run= opt[["dry_run"]])
} 
if( FALSE)
{
#### Step 2: Create the man_review.txt file by reviewing the results, following http:// ....
   set_allelic_funcs()
   apply_review_and_extract( pp.review.fn=calls_FN, obj.name=obj.name, analyst.id="SLC" )
}

