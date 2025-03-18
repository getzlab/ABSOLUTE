## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

ApplySSNVModel <- function(mode.res, mut.cn.dat, SSNV_model, verbose=FALSE)
{
  if (verbose) {
    print(paste("Evaluating ", nrow(mut.cn.dat),
                " mutations over ", nrow(mode.res[["mode.tab"]]),
                " purity/ploidy modes: ", 
                sep = ""))
  }

  n_modes = nrow(mode.res[["mode.tab"]])
  ccf_grid = SSNV_model[["ccf_grid"]]
  mode.res[["SSNV.ccf.dens"]] = array( NA, dim=c(n_modes, nrow(mut.cn.dat), length(ccf_grid)) )
  dimnames(mode.res[["SSNV.ccf.dens"]])[[3]] = ccf_grid
  
  mode.res[["modeled.muts"]] = list()
  mode.res[["mode_SSNV_models"]] = list()

  for (j in 1:n_modes) 
  {
    alpha <- mode.res[["mode.tab"]][j, "alpha"]
    seg.q.tab <- mode.res[["mode_SCNA_models"]][[j]] [["seg.q.tab"]]
    
    subclonal_scna_tab = mode.res[["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]][j, , ]
    SCNA_log_ccf_dens = mode.res[["subclonal_SCNA_res"]][["log_CCF_dens"]][j, , ]
    
    res = FitPPModeSomaticMuts(SSNV_model, mut.cn.dat, mode.res$mode.tab[j, ], subclonal_scna_tab, SCNA_log_ccf_dens, seg.q.tab, verbose=verbose)
        
    modeled.muts <- res[["modeled.muts"]]
    som.theta.q.map <- res[["som.theta.q.map"]]
    mode.res[["SSNV.ccf.dens"]][j,,] <- res[["ccf.dens"]] 
    mode.res[["modeled.muts"]][[j]] <- cbind(modeled.muts, purity = alpha, SSNV_skew=SSNV_model[["SSNV_skew"]] )

    mode.res[["mode_SSNV_models"]][[j]] = res[["mode_SSNV_models"]]
    
    mode.res[["mode.tab"]][j, "SSNV_LL"] <- sum(modeled.muts[, "LL"], na.rm = TRUE) + dDirichlet(som.theta.q.map, SSNV_model[["kPiSomThetaQ"]], log.p = TRUE)
 
    if (verbose) { cat(".") }
  }
  if (verbose) {
    cat("\n")
  }
  
#  new.ll <- mode.res[["mode.tab"]][, "combined_LL"] +
#            mode.res[["mode.tab"]][, "SSNV_LL"]
#  mode.res[["mode.tab"]][, "combined_LL"] <- new.ll
  
  return(mode.res)
}




FitPPModeSomaticMuts <- function(SSNV_model, mut.cn.dat, mode_info, subclonal_scna_tab, 
			   scna_log_ccf_dens, seg.q.tab, verbose=FALSE) 
{
  clonal_scna_mut_ix = !get_subclonal_scna_mut_ix(mut.cn.dat, subclonal_scna_tab)

  res = get_subclonal_scna_tab(mut.cn.dat[!clonal_scna_mut_ix,], subclonal_scna_tab, seg.q.tab )
  subclonal_scna_tab = res[["subclonal_scna_tab"]]
  subclonal.mut.tab = as.data.frame(matrix( NA, nrow=nrow(mut.cn.dat), ncol=ncol(res[["subclonal.mut.tab"]]) ))
  colnames(subclonal.mut.tab) = colnames(res[["subclonal.mut.tab"]])

  subclonal.mut.tab[!clonal_scna_mut_ix,] = res[["subclonal.mut.tab"]]
  clonal.mut.tab = get_muts_nearest_clonal_scna(mut.cn.dat, seg.q.tab, SSNV_model[["kQ"]])
  mut.modeled.cn = cbind(clonal.mut.tab, subclonal.mut.tab, "clonal_scna_mut_ix"=clonal_scna_mut_ix)

  fit_res = fit_SSNV_model( cbind(mut.cn.dat, mut.modeled.cn), mode_info, SSNV_model, subclonal_scna_tab, scna_log_ccf_dens )
  som_theta_q_map = fit_res$som_theta_Q_MAP
  post_prs = fit_res$post_Prs
  ssnv.ccf.dens = fit_res[["ssnv.ccf.dens"]]
     
  ## Subclonal SCNA?
  var_classes = ClassifySomaticVariants(post_prs, 0.5)
  q_s = post_prs[, "modal_q_s"]
  mut_mult_res = calc_CCF_95CI( cbind(mut.cn.dat, mut.modeled.cn), ssnv.ccf.dens, mode_info, q_s, post_prs[, "Pr_somatic_clonal"], SSNV_model )
  ##

  detection_power = mode_SSNV_pow_calc( SSNV_model, mut.cn.dat, mut.modeled.cn, mode_info["alpha"] )
  detection_power_for_single_read = mode_SSNV_pow_calc( SSNV_model, mut.cn.dat, mut.modeled.cn, mode_info["alpha"], single_read=TRUE )

  modeled.muts <- cbind( mut.modeled.cn, post_prs, var_classes, mut_mult_res, fit_res[["H.ev"]], detection_power=detection_power, detection_power_for_single_read=detection_power_for_single_read )
#  modeled.muts <- cbind( mut.modeled.cn, post_prs, fit_res[["H.ev"]] )

  return(list( "ccf.dens" = ssnv.ccf.dens, "modeled.muts" = modeled.muts, "som.theta.q.map" = som_theta_q_map, "mode_SSNV_models"=fit_res[["SSNV_model"]]))
}






CreateMutCnDat <- function(maf, indel.maf, seg.dat, min.mut.af, verbose=FALSE) 
{
  indel_filters = function(maf)
  {
     class = maf[, "Variant_Classification"]
     n.ix = class %in% c("IGR", "Intron", "5'UTR", "3'UTR", "RNA", "5'Flank")

     msg = paste("Removing ", sum(n.ix), " of ", length(n.ix), " indels in IGR / Intron / UTR / Flank / RNA", sep="")
     print(msg)
     maf = maf[!n.ix,]

## Turn off "don't ask don't tell" for indels; seems hard to hallucinate these reads due to seq errors
     if( FALSE & nrow(maf) > 0)
     {
#        ix = is.na(maf[,"i_judgement"])  ## only here due to forced-calling
#        ix = maf[,"i_judgement"] == "REJECT"  ## only here due to forced-calling
        ix = maf[,"force_called_site"] == "YES"
        ix[is.na(ix)] = FALSE 
        maf[ix,"t_alt_count"] = 0
     }

#     type = maf[, "Variant_Type"] ## "INS" or "DEL"
     return(maf)
  }


  if( !is.null(indel.maf) )
  {
     indel_missing_cols = setdiff(colnames(maf), colnames(indel.maf))
     imc = matrix( NA, nrow=nrow(indel.maf), ncol=length(indel_missing_cols))
     colnames(imc)= indel_missing_cols
     indel.maf = cbind(indel.maf, imc)
     nc = intersect(colnames(indel.maf), colnames(maf)) 
     filtered.indel.maf = indel_filters(indel.maf)
     maf = rbind(maf, filtered.indel.maf[,nc] )
  }

  mut.cn.dat <- maf
  
  if ("total_normals_called" %in% colnames(mut.cn.dat)) {
    ix <- mut.cn.dat[, "total_normals_called"] > 1
    if (verbose) {
      print(paste("Removing ", sum(ix), " of ", length(ix),
                  " mutations due to seen in > 1 normals", 
                  sep = ""))
    }
    mut.cn.dat <- mut.cn.dat[!ix, ]
  }
  
  if ("dbSNP_Val_Status" %in% colnames(mut.cn.dat)) {
    mut.cn.dat[["dbSNP_Val_Status"]][is.na(mut.cn.dat[["dbSNP_Val_Status"]])] <- ""
  }
  
  cols <- colnames(mut.cn.dat)
  
  cix <- which(cols %in% c("i_t_ref_count", "t_ref_count"))
  colnames(mut.cn.dat)[cix] <- c("ref")
  
  cix <- which(cols %in% c("i_t_alt_count", "t_alt_count"))
  colnames(mut.cn.dat)[cix] <- c("alt")
  
  cix <- which(cols %in% c("dbSNP_Val_Status"))
  colnames(mut.cn.dat)[cix] <- "dbSNP"
  
  cix <- which(cols %in% c("Tumor_Sample_Barcode"))
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "sample"
  }

  if( "failure_reasons" %in% cols ) {
     colnames(mut.cn.dat)[ which(colnames(mut.cn.dat)=="failure_reasons") ] = "i_failure_reasons" 
     cols[ which(cols=="failure_reasons") ] = "i_failure_reasons" 
  }

## save actual observed count of alt reads
  mut.cn.dat[["observed_alt"]] = mut.cn.dat[["alt"]]

## Turn on forced-calling of SSNVs
  if( FALSE & "i_failure_reasons" %in% cols )
  {
    ix1 = grep( "fstar_tumor_lod", mut.cn.dat[, "i_failure_reasons"] )
    ix2 = which(mut.cn.dat[,"alt"] > 0)
    ix = intersect( ix1, ix2)
    
    if( length(ix) > 0 ) {
       mut.cn.dat[ix, "alt"] = 0
    }

    if( verbose ) {
       msg = paste( length(ix), " mutations rejected for fstar_tumor_lod, alt set to 0", sep="" )
       print(msg)
    }

    is.important = identify_potential_clinically_actionable_mutations( mut.cn.dat )
    save.ix = is.important & mut.cn.dat[,"observed_alt"] > mut.cn.dat[,"alt"] 
 
    if( any(save.ix) )
    {
       mut.cn.dat[save.ix,"alt"] = mut.cn.dat[save.ix, "observed_alt"]
    }
    if(verbose) 
    {
       print( paste( "Reverting ", sum(save.ix), " important mutations to force-called alt reads:", sep=""))
       if( any(save.ix) ) { print(mut.cn.dat[save.ix, c("Hugo_Symbol", "Variant_Classification")] )  }
    }
  }

  
#  na.ix <- apply(is.na(mut.cn.dat[, c("ref", "alt")]), 1, sum) > 0
  na.ix <- is.na(mut.cn.dat[, "alt"] + mut.cn.dat[, "ref"])
  if (verbose) {
    print(paste("Removing ", sum(na.ix), " of ", length(na.ix),
                " mutations with NA coverage", 
                sep = ""))
  }
  mut.cn.dat <- mut.cn.dat[!na.ix, ]

## check for negative read-counts (yes, this happens sometimes for indels)
  neg.ix = mut.cn.dat[, "alt"] < 0 
  if (verbose) { print(paste("Setting ", sum(neg.ix), " of ", length(neg.ix), " alt read counts < 0 to 0", sep = "")) }
  mut.cn.dat[neg.ix,"alt"] = 0
  neg.ix = mut.cn.dat[, "ref"] < 0
  if (verbose) { print(paste("Setting ", sum(neg.ix), " of ", length(neg.ix), " ref read counts < 0 to 0", sep = "")) }
  mut.cn.dat[neg.ix,"ref"] = 0
  
  af <- mut.cn.dat[, "alt"] / (mut.cn.dat[, "alt"] + mut.cn.dat[, "ref"])
  ix <- af < min.mut.af
  ix[  (mut.cn.dat[, "alt"] + mut.cn.dat[, "ref"]) == 0 ] = FALSE 

  if (verbose) {
    print(paste("Removing ", sum(ix), " of ", length(ix),
                " mutations due to allelic fraction < ", 
                min.mut.af, sep = ""))
  }

  if (sum(!ix) == 0) {  
    stop("no mutations left!") 
  }
  mut.cn.dat = mut.cn.dat[!ix, , drop=FALSE]
    
  mut.seg.ix <- GetMutSegIx(mut.cn.dat, seg.dat[["segtab"]])  
  ix <- apply(is.na(mut.seg.ix), 1, sum) == 0

## support HSCN
  male_X = seg.dat[["obs.scna"]][["male_X"]][mut.seg.ix[,1]]
  male_X[ !ix ] = FALSE

  mut.cn.dat = data.frame( mut.cn.dat, "male_X"=male_X )

  if( any(is.na(male_X))){ stop() }

  if (verbose) {
    print( paste( "Mapped ", sum(male_X), " mutations with male_X status", sep=""))
    print(paste("Removing ", sum(!ix), " unmapped mutations on Chrs: ", sep = ""))
    print(mut.cn.dat[!ix, "Chromosome"])
  }
  if (sum(ix) == 0) {
    stop("No mutations left")
  }

  mut.cn.dat <- mut.cn.dat[ix, ]
  mut.seg.ix <- mut.seg.ix[ix, , drop = FALSE]
  mut.cn.dat <- cbind(mut.cn.dat, mut.seg.ix)

  mut.cn.dat = select_protein_change_annot_using_COSMIC( mut.cn.dat, verbose=verbose )
  
  return(mut.cn.dat)
}




identify_potential_clinically_actionable_mutations = function( maf )
{
   #data("VanAllen2014_TARGET", package="ABSOLUTE")  ## provides TARGET
   target_genes = TARGET

   silent_classes = c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region", "De_novo_Start_InFrame") 

#   LOF_classes = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation") 
   crit_col = "Types_of_recurrent_alterations"

   mut_genes.ix = grep( "Mutation", target_genes[, crit_col] )
   biallelic_genes.ix = grep( "Biallelic Inactivation", target_genes[,crit_col] )
   gene.ix = union(mut_genes.ix, biallelic_genes.ix)
   genes= target_genes[gene.ix, "Gene"]
   
   is.important = maf[,"Hugo_Symbol"] %in% genes & !(maf[,"Variant_Classification"] %in% silent_classes)

   return(is.important)
}


## Selects between UNIPROT or other AA change using COSMIC
select_protein_change_annot_using_COSMIC = function( MAF, verbose=FALSE )
{
   get_COSMIC_count = function( Hugo_Symbol, Protein_Change, pankey_counts )
   {
# count # of identical codon changes in COSMIC
      keys = paste(Hugo_Symbol, "__", Protein_Change, sep="")
      var_count = rep(0, length(keys) )
      ix = which( keys %in% names(pankey_counts) )
      var_count[ix] = pankey_counts[ keys[ix] ]

      return(var_count)
   }

   #data("COSMIC_protein_change_counts_v67_241013", package="ABSOLUTE")
   #load("/xchip/tcga/Tools/absolute/releases/v1.5/data/COSMIC_protein_change_counts_v67_241013.RData")
   pankey_counts = COSMIC_protein_change_counts

# count # of identical codon changes in COSMIC
   keys = paste( MAF[, "Hugo_Symbol"], "__", MAF[, "Protein_Change"], sep="")
   other_count = get_COSMIC_count( MAF[, "Hugo_Symbol"], MAF[, "Protein_Change"], pankey_counts )


   pc = MAF[,"Protein_Change"]
   pc[is.na(pc)] = ""
   MAF[,"Protein_Change"] <- as.character(pc)

   res = strsplit( MAF[,"Protein_Change"], "[0-9]+" )
   A1 = unlist( lapply( res, "[", 1 ))
   A2 = unlist( lapply( res, "[", 2 ))
   UNIPROT_Protein_Change = paste( A1, MAF[,"UniProt_AApos"], A2, sep="" )

   uniprot_count = get_COSMIC_count( MAF[, "Hugo_Symbol"], UNIPROT_Protein_Change, pankey_counts )

   switch_to_uniprot =  uniprot_count >= other_count & MAF[,"Protein_Change"] != "" & !is.na(MAF[,"UniProt_AApos"]) #always uniprot

   if( verbose ) 
   {
      print( paste("Switching to UNIPROT for ", sum(switch_to_uniprot), " mutations based on COSMIC counts ", sep=""))

     
      print( cbind( MAF[ switch_to_uniprot, "Hugo_Symbol"], 
                    MAF[ switch_to_uniprot, "Protein_Change"],
                    other_count[ switch_to_uniprot ] ) )

     
      print( cbind( MAF[ switch_to_uniprot, "Hugo_Symbol"], 
                    UNIPROT_Protein_Change[switch_to_uniprot],
                    uniprot_count[ switch_to_uniprot ] ))
   }
   MAF[ switch_to_uniprot, "Protein_Change"] = UNIPROT_Protein_Change[switch_to_uniprot]
   COSMIC_count = other_count
   COSMIC_count[ switch_to_uniprot ] = uniprot_count

   MAF = cbind( MAF, "Number of times codon change is in COSMIC"=COSMIC_count )

   return(MAF)
}

