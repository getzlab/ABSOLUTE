source("~scarter/CGA/R/broad_utils/LSF_scatter.R")


#batch_exec_ABSOLUTE = function( ABSOLUTE_argv, obj.name, MAF_DIR, SIF_FN, sample_list_FN, MAF_SIF_FN, queue="hour", overwrite=FALSE, dry_run=FALSE )
batch_exec_ABSOLUTE = function( ABSOLUTE_argv, obj.name, var_bsub_argv, queue="hour", overwrite=FALSE, dry_run=FALSE, num_solutions_plotted=3 )
{
   ABSOLUTE_argv = insert_default_args( ABSOLUTE_argv )

  ## setup LSF exec control args
   control_argv = list()
	control_argv$OVERWRITE= overwrite
	control_argv$DRY_RUN = dry_run
	control_argv$wait = TRUE
	control_argv$QUEUE = queue
	control_argv$R_STUB_FN = file.path(CGA_DIR, "broad_utils/ABSOLUTE_stub.R" )
        control_argv$BJOB = obj.name

        OUT_DIR_base = "." #file.path( "UTE_results", obj.name )
        control_argv$R_DIR = file.path( OUT_DIR_base, "R/" )

        RESULTS_DIR =  file.path( OUT_DIR_base  )
        ABSOLUTE_argv$results.dir = file.path( "results" )

   print ("try lsf")
   #LSF_scatter( control_argv, ABSOLUTE_argv, var_bsub_argv )


 #### Step 2: collect results for summarization:
   file.base = paste( var_bsub_argv[, "output.fn.base"], ".ABSOLUTE", sep = "")
   absolute.files= file.path( ABSOLUTE_argv$results.dir, paste(file.base, "RData", sep = "."))

print (absolute.files)

   CreateReviewObject( obj.name, absolute.files, RESULTS_DIR, ABSOLUTE_argv$copy_num_type, plot.modes=TRUE, num_solutions_plotted=num_solutions_plotted, verbose=TRUE) 
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



insert_default_args = function( ABSOLUTE_argv )
{
   default_abs_args = list( min.ploidy=1, max.ploidy=8,
                        max.as.seg.count=1500, max.non.clonal=0.5, max.neg.genome=0.05,
                        maf.fn=NA, indel.maf.fn=NA, min.mut.af = 0,
                        output.fn.base=NA, min_probes=10, max_sd=100, sigma.h=0.01, 
			SSNV_skew=1, filter_segs=TRUE, force.alpha=NA, force.tau=NA, verbose=TRUE )

   for( i in 1:length(default_abs_args) )
   {
      if( !( names(default_abs_args)[i] %in% names(ABSOLUTE_argv) ) ) 
      {
         ABSOLUTE_argv[[  names(default_abs_args)[i]  ]] = default_abs_args[[i]] 
      }
   }

   return( ABSOLUTE_argv )
}

