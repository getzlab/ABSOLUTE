## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


## FIXME: Some funcs don't seem to be used at all:
##   - GetMutsScna
##   - CalcSampleMutsPostPr
##   - get_cr_grid_from_ccf_grid

## TODO:  Some of these are for general allelic analysis.  Others are for specific HAPSEG / SNP error models.  Sort out which are which.

set_allelic_funcs = function() {

## Allelic models
  ExtractSampleObs <<- AllelicExtractSampleObs

  GetMutSegIx <<- AllelicGetMutSegIx
  get_muts_nearest_clonal_scna <<- allelic_get_muts_nearest_clonal_scna
  CalcClonalSSNVLoglik <<- AllelicCalcClonalSSNVLoglik
  CalcGermlineSNVLoglik <<- AllelicCalcGermlineSNVLoglik

  CalcChrArmDistr <<- AllelicCalcChrArmDistr
  GetCopyRatioComb <<- AllelicGetCopyRatioComb
  get_cr_grid_from_ccf_grid <<- allelic_get_cr_grid_from_ccf_grid
  GetAbsSegDat <<- AllelicGetAbsSegDat
  calc_sample_muts_on_subclonal_scna <<- allelic_calc_sample_muts_on_subclonal_scna

  get_subclonal_scna_mut_ix <<- allelic_get_subclonal_scna_mut_ix  
  get_subclonal_SCNA_info  <<- allelic_get_subclonal_SCNA_info
  get_subclonal_scna_tab <<- get_allelic_subclonal_scna_tab
  annotate_SSNVs_on_clonal_SCNAs <<- allelic_annotate_SSNVs_on_clonal_SCNAs
  eval_SNV_models_evidence <<- allelic_eval_SNV_models_evidence

## HAPSEG SNP model
  get_seg_sigma <<- Allelic_get_seg_sigma
} 

set_total_funcs = function() {

  ExtractSampleObs <<- total_extract_sample_obs
  GetMutSegIx <<- total_get_mut_seg_ix    
  get_muts_nearest_clonal_scna<<- total_get_muts_nearest_clonal_scna
  CalcChrArmDistr <<- total_calc_chr_arm_distr
  GetCopyRatioComb <<- total_get_copy_ratio_comb
  get_cr_grid_from_ccf_grid <<- total_get_cr_grid_from_ccf_grid
  GetAbsSegDat <<- total_get_abs_seg_dat

#  calc_sample_muts_on_subclonal_scna <<- total_calc_sample_muts_on_subclonal_scna  
  get_subclonal_scna_mut_ix <<- total_get_subclonal_scna_mut_ix
  get_subclonal_SCNA_info <<- total_get_subclonal_SCNA_info

  get_subclonal_scna_mut_ix <<- total_get_subclonal_scna_mut_ix 
  get_subclonal_scna_tab <<- get_total_subclonal_scna_tab
  annotate_SSNVs_on_clonal_SCNAs <<- total_annotate_SSNVs_on_clonal_SCNAs
  eval_SNV_models_evidence <<- total_eval_SNV_models_evidence

  get_seg_sigma <<- CAPSEG_get_seg_sigma
}
