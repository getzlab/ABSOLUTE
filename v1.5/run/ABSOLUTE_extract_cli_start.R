library(optparse)
option.list <- list(
  make_option("--solution_num", dest="solution_num"),
  make_option("--analyst_id", dest="analyst_id"),
  make_option("--rdata_modes_fn", dest="rdata_modes"),
  make_option("--sample_name", dest="sample_name"),
  make_option("--results_dir", dest="results_dir",
              default=getwd()),
  make_option("--abs_lib_dir", dest="abs_lib_dir"),
  make_option("--max_ploidy", dest="max_ploidy",type="double",default=6.0))

#--seg_dat_fn --maf_fn --indelmaf_fn --sample_name --results_dir --ssnv_skew

#"""Rscript /run/ABSOLUTE_cli_start.R \
# --seg_dat_fn {seg_file} --maf_fn {snv_maf} --indelmaf_fn {indel_maf} --sample_name {sample_name} --results_dir results \
#--ssnv_skew {skew} """

# --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

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
#indel.maf.fn= opt[["indelmaf_fn"]]
solution_num= opt[["solution_num"]]
min.mut.af= 0
output.fn.base= opt[["sample_name"]]
max_sd= 100
modesegs.fn= opt[["rdata_modes"]]
filter_segs= TRUE
force.alpha= NA
force.tau= NA
verbose= TRUE
results.dir= opt[["results_dir"]]
sample.name= opt[["sample_name"]]
analyst.id= opt[["analyst_id"]]
gender= NA
#library(ABSOLUTE)



require(GenomicRanges)
require(gplots)
require(RColorBrewer)

#suppressPackageStartupMessages(require(gplots))
#suppressPackageStartupMessages(require(RColorBrewer))
CGA_DIR_ABS=opt[["abs_lib_dir"]] ##"/soft/local/absolute" #/xchip/tcga/Tools/absolute/releases/v1.5/

   print( paste("sourcing files in ", CGA_DIR_ABS, sep=""))
   #rr = dir( file.path( CGA_DIR_ABS, "ABSOLUTE/sandbox"), full.names=TRUE )
   rr = dir(CGA_DIR_ABS,full.names=TRUE, pattern = '.R$')
   for( i in 1:length(rr) ) { 
   source(rr[i])}


# moved into the absolute root dir
#   source( file.path( CGA_DIR_ABS, "Phylogic", "log_DP_post_fixed.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "log_DP_post.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "run_DP.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "DP_utils.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "summarize_DP.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "CCF_DP_plots.R" ) )



#debug (RunAbsolute)
load(file.path( CGA_DIR_ABS,"data","ChrArmsDat.RData"))
load(file.path( CGA_DIR_ABS,"data","ChrArmPriorDb.RData"))
load(file.path( CGA_DIR_ABS,"data","diseaseMap.RData"))
#RunAbsolute( seg.dat.fn, primary.disease, platform, sample.name, results.dir, copy_num_type,
#                        genome_build, gender, min.ploidy, max.ploidy,
#                        max.as.seg.count, max.non.clonal, max.neg.genome,
#                        maf.fn, indel.maf.fn, min.mut.af,
#                        output.fn.base, min_probes, max_sd, sigma.h, SSNV_skew,
#			filter_segs, force.alpha, force.tau, verbose )


file.base = paste( output.fn.base, ".ABSOLUTE", sep = "")
absolute.files= file.path( results.dir, paste(file.base, "RData", sep = "."))
print (absolute.files)
#CreateReviewObject( sample.name, absolute.files, ".", copy_num_type, plot.modes=TRUE, num_solutions_plotted=num_solutions_plotted, verbose=TRUE)

set_allelic_funcs()
called.segobj.list = run_PP_calls_liftover_from_num(as.numeric(solution_num), analyst.id, modesegs.fn, results.dir, sample.name, verbose=verbose )

# if object found - OK, else fail
if( length(called.segobj.list) > 0 )
   {
      ExtractReviewedResults( called.segobj.list, analyst.id , results.dir, sample.name, verbose=TRUE )
   } else { 
stop("called.segobj.list has length 0!") 
}

