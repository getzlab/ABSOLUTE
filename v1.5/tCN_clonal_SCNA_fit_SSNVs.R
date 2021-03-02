## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


total_get_mut_seg_ix = function(maf, segtab) {
  ## compute lookup table for each mutation into seg_Q_tab
  N = nrow(maf)
  mut_seg_ix = matrix(NA, nrow=N, ncol=1)
  colnames(mut_seg_ix) = c("mut_seg_ix")
  
  ##  This is Dan-Avi's fault.
  if (!("End_position" %in% colnames(maf))) {
    end_position = maf[, "Start_position"]
    maf = cbind(maf, End_position=end_position)
  }
  
  for (i in seq_len(N)) {
    seg.ix = maf[i, "Chromosome"] == segtab[, "Chromosome"] &
      maf[i, "Start_position"] >= segtab[, "Start.bp"] &
      maf[i, "End_position"] <= segtab[, "End.bp"] 
    
    if (sum(seg.ix) != 1) { 
      next 
    }
    
    seg_id = which(seg.ix)
    mut_seg_ix[i, 1] = seg_id
  }
## TODO: map missing SSNVs to closest adjacent seg
print("TODO: map missing SSNVs to closest adjacent seg!")
  
#  if( any(is.na(mut_seg_ix)) ) { stop() }
  
  return(mut_seg_ix)
}

get_tcn_somatic_mut_comb = function(alpha, q) {
  q_a = c(1:q)
  
  f.reads = alpha * q_a / (2 * (1 - alpha) + alpha * q)
  vals = q_a
  res = cbind(vals, f.reads)
  
  return(res)
}

get_tcn_somatic_X_mut_comb = function(alpha, q) {
  q_a = c(1:q)
  
  f.reads = alpha * q_a / (1 * (1 - alpha) + alpha * q)
  vals = q_a
  res = cbind(vals, f.reads)
  
  return(res)
}


get_tcn_germline_mut_comb = function(alpha, q) {
  eps = 1e-3
  q_vals = c(0:q)
  
  f.reads = ((1 - alpha) + alpha * q_vals) / (2 * (1 - alpha) + alpha * q)   
  f.reads[3] = 1 - eps
  vals = q_vals
  res = cbind(vals, f.reads)
  
  return(res)
}

get_tcn_germline_X_mut_comb = function(alpha, q) {
  eps = 1e-3
  q_vals = c(0:q)
  
  f.reads = (alpha * q_vals) / ((1 - alpha) + alpha * q)   
  f.reads[3] = 1 - eps
  vals = q_vals
  res = cbind(vals, f.reads)
  
  return(res)
}


TotalCalcClonalSSNVLoglik <- function(mut.cn.dat, mode_info, SSNV_model)
{
  som.theta.q = SSNV_model[["som_theta_Q_mode"]]
  mut_class_w = SSNV_model[["mut_class_w"]]

  LL_mat <- matrix(-Inf, nrow=nrow(mut.cn.dat), ncol=length(som.theta.q))

  alpha = mode_info["alpha"]

  u.keys = unique(mut.cn.dat[, "q_hat"])
  n.keys = length(u.keys)
  mut_n.keys = rep(NA, n.keys)
  
  for (i in seq_len(n.keys)) {
    qt = u.keys[i]
    mut.ix = mut.cn.dat[, "q_hat"] == qt 
    mut_n.keys[i] = sum(mut.ix)
  }
  u.keys = u.keys[order(mut_n.keys, decreasing=TRUE)]
  mut_n.keys = sort(mut_n.keys, decreasing=TRUE)
 

  for (i in seq_len(n.keys)) 
  {
    qt = u.keys[i]
    mut.ix = mut.cn.dat[, "q_hat"] == qt 

    if (qt==0) {
      LL_mat[mut.ix, ] <- -Inf
      next
    }

    alt = mut.cn.dat[mut.ix, "alt"]
    ref =  mut.cn.dat[mut.ix, "ref"]

## TODO: add support for X-chr
#    m.ix = mut.cn.dat[mut.ix, "male_X"]

    ## Somatic autosome model calc
    som.comb  = get_tcn_somatic_mut_comb( alpha, qt )
    som.w <- clonal_SSNV_Mult_Prior(mut_class_w, nrow(som.comb), som.theta.q)
    som.w <-  matrix(som.w, ncol=nrow(som.comb), nrow=length(alt), byrow=TRUE)
    som.ll <- log(som.w) + SSNV_FhatCombPost( alt, ref, som.comb[,"f.reads"], SSNV_model )

    LL_mat[mut.ix, c(1:ncol(som.ll))] = som.ll
  }

  return(LL_mat)
}   


# version for tCR clonal SCNA segs
## data are all from a single sample
total_eval_SNV_models_evidence = function(mut.cn.dat, mode_info, post.ccf.LL.grid, SSNV_model ) 
{
  ## n_mut X Q 
   clonal.comb.LL = TotalCalcClonalSSNVLoglik(mut.cn.dat, mode_info, SSNV_model)
   Z = LogAdd(clonal.comb.LL)
   Z[is.nan(Z)] = -Inf
   som_mut_Q_tab = exp(clonal.comb.LL - Z)
   homdel.ix = mut.cn.dat[, "q_hat"]==0

  ## n_mut X 1
   som.clonal.log.ev = LogAdd(clonal.comb.LL) #+ log(SSNV_model[["mut_class_w"]][["SM"]])
   som.clonal.log.ev[is.nan(som.clonal.log.ev)] = -Inf

   subclonal.exp.log.ev = H_exp_subclonal_SSNV( post.ccf.LL.grid, SSNV_model )
   subclonal.unif.log.ev = H_unif_subclonal_SSNV( post.ccf.LL.grid, SSNV_model )

   if( SSNV_model[["mut_class_w"]][["GL"]] > 0 )
   {
      stop( "Germline multiplcity not implemented for tCR input" )
   } else {
      germline.log.ev= -Inf
   }

   SSNV.on.CN0.log.ev = rep(-Inf, nrow(mut.cn.dat))
   SSNV.on.CN0.log.ev[ homdel.ix ] = SSNV_model[["SSNV.on.CN0"]]

## put together mutually exclusive hypothoses:
   mut.log.ev.mat = cbind( som.clonal.log.ev, subclonal.exp.log.ev, subclonal.unif.log.ev, germline.log.ev, SSNV.on.CN0.log.ev )

   mut.log.ev.mat[ homdel.ix, c(1,2,3) ] = -Inf  ## for hom dels
   LL = LogAdd(mut.log.ev.mat)
   mut.ev.mat = exp(mut.log.ev.mat - LL)

   post_Prs = total_annotate_SSNVs_on_clonal_SCNAs( mut.ev.mat, som_mut_Q_tab, mut.cn.dat, SSNV_model, LL)

   return(list("mut.ev.mat"=mut.ev.mat, "som_mut_Q_tab"=som_mut_Q_tab, "post_Prs"=post_Prs))
}


total_annotate_SSNVs_on_clonal_SCNAs = function( mut.ev.mat, som_mut_Q_tab, mut.cn.dat, SSNV_model, LL )
{
   mut.w = SSNV_model[["mut_class_w"]]

# Original contract for CalcSampleMutsPostPr()
   cols <- c("Pr_somatic_clonal", "Pr_germline", "Pr_subclonal",
             "Pr_subclonal_wt0", "Pr_wt0", "Pr_ge2", "Pr_GL_som_HZ_alt",
             "Pr_GL_som_HZ_ref", "Pr_cryptic_SCNA", "modal_q_s", "LL")
   post_Prs = matrix( 0, nrow=nrow(mut.cn.dat), ncol=length(cols) )
   colnames(post_Prs)=cols
   post_Prs[,"LL"] = LL

   LOH.ix = mut.cn.dat[, "q_hat"]==1  
   if( any(LOH.ix) & mut.w[["GL"]] > 0) 
   {
      GL_PR = rowSums(mut.ev.mat[LOH.ix, c("alt_minor", "alt_major", "hom_alt"), drop=FALSE])
      post_Prs[LOH.ix, "Pr_GL_som_HZ_alt"] <- mut.ev.mat[LOH.ix, "alt_major"] / GL_PR
      post_Prs[LOH.ix, "Pr_GL_som_HZ_ref"] <- mut.ev.mat[LOH.ix, "alt_minor"] / GL_PR
   }

   post_Prs[,"Pr_somatic_clonal"] = mut.ev.mat[,1]
   post_Prs[,"Pr_ge2"] = mut.ev.mat[,1] * rowSums(som_mut_Q_tab[,-1, drop=FALSE])

   post_Prs[,"Pr_subclonal"] = rowSums(mut.ev.mat[,c(2,3), drop=FALSE])
   post_Prs[,"Pr_germline"] = mut.ev.mat[,4]

   if( any(LOH.ix))
   {
      s1 = som_mut_Q_tab[LOH.ix,, drop=FALSE ]
      r = (1:nrow(s1)-1)
      q2 = mut.cn.dat[LOH.ix, "q_hat"]
      post_Prs[LOH.ix,"Pr_wt0"] = t(s1)[ q2 + ncol(s1) * r ]  * post_Prs[LOH.ix, "Pr_somatic_clonal"] 
   }

   ix.01 = LOH.ix
   post_Prs[ix.01, "Pr_subclonal_wt0"] = post_Prs[ix.01, "Pr_subclonal"]

#if( any(ix.01 & post_Prs[,"Pr_wt0"] != post_Prs[,"Pr_somatic_clonal"]) ) { stop() }

  modal_q_s <- rep(NA, nrow(som_mut_Q_tab))

 ## assume subclonal muts have mult 1 in SC fraction.
  modal_q_s <- rep(NA, nrow(som_mut_Q_tab))
  if (mut.w[["SM"]] > 0) {
    nix <- post_Prs[, "Pr_somatic_clonal"] < 0.5
    nix[is.na(nix)] = TRUE
    if (any(!nix)) {
      modal_q_s[!nix] <- apply(som_mut_Q_tab[!nix, , drop=FALSE], 1, which.max)
    }
    modal_q_s[nix] <- 1
  }
  post_Prs[,"modal_q_s"]=modal_q_s


  if(any(is.na( post_Prs[,"Pr_somatic_clonal"] ))) { stop("NA in Pr_somatic_clonal") }

  return(post_Prs)
}



