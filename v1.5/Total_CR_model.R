## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

total_get_copy_ratio_comb = function(Q, delta, b, error_model) {
  xx = (delta * (c(1:Q) - 1) + b)
  return(xx)
}


total_get_cr_grid_from_ccf_grid = function(qc, qs, delta, comb, ccf_grid) {
  d = qs - qc
  d_d = delta * d
  cr_grid = (d_d * ccf_grid) + comb[qc + 1]
  
  if (d < 0) { 
    cr_grid = rev(cr_grid)
  }
  
  return(cr_grid)
}


get_tcr_chr_arm_segs = function(seg.obj, chr_arm_dat) {
  seg.dat = seg.obj$segtab
  
  int_w = c()
  ix = c()
  
  arm_len_bp = chr_arm_dat["End.bp"] - chr_arm_dat["Start.bp"]
  
  for (i in 1:nrow(seg.dat)) {
    if( seg.dat[i,"Chromosome"] != chr_arm_dat["chr"]) { 
      next 
    }
    int_start = max(as.numeric(c(seg.dat[i, "Start.bp"], chr_arm_dat["Start.bp"])))
    int_end = min(as.numeric(c(seg.dat[i, "End.bp"], chr_arm_dat["End.bp"])))
    
    ## seg does not overlap region
    if (int_start > int_end) { 
      next 
    }
    int_len_bp = int_end - int_start
    int_w = c(int_w, int_len_bp / arm_len_bp)
    
    ix = c(ix, i)
  }
  
  int_w = as.numeric(int_w)
  
  return(list("int_W"=int_w, "ix"=ix ))
}


total_calc_chr_arm_distr = function(seg.obj, seg_q, chr_arms_dat) {
  n_arm = nrow(chr_arms_dat)
  chr_arm_tab = array(NA, dim=c(1, n_arm, ncol(seg_q)))
  
  for (i in seq_len(n_arm)) {
    chr_dat = get_tcr_chr_arm_segs(seg.obj, chr_arms_dat[i, ])      
    
    if (length(chr_dat$int_W) == 0 ) { 
      next 
    }
    
    chr_arm = array(NA, dim=c(1, ncol(seg_q)))
    chr_arm[1, ] = colSums(seg_q[chr_dat$ix, , drop=FALSE] * chr_dat$int_W)
    
    chr_arm_tab[, i, ] = 0
    chr_arm_tab[1, i, which.max(chr_arm[1,])] = 1
  }
  
  return(chr_arm_tab)
}


## This version is for data from raw CAPSEG - no seg-level std errs!
total_get_abs_seg_dat = function(segobj) {
  Q = 15
  qq = Q
  
  # Get column number of the max of each row and the expected 
## ONLY USE THESE FIELDS FROM MODEL: 
#  seg_qz_tab = segobj$mode.res$seg.qz.tab[1, , ]
#  seg_q_tab = segobj$mode.res$seg.q.tab[1, , ]

  seg_q_tab <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["seg.q.tab"]]
  seg_qz_tab <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["seg.qz.tab"]]

##
  
  max_mat <- apply(seg_qz_tab, MARGIN=1, function(x) which.max(x))
  subclonal_ix = (max_mat == (Q + 1))
  
  max_mat = apply(seg_q_tab, MARGIN=1, function(x) which.max(x))
  
  exp_mat = apply(seg_q_tab, MARGIN=1, 
                  function(x) {
                    x <- x[1:qq] / sum(x[1:qq])
                    return(sum(x * c(1:qq))) 
                  })


  
  # seg_list is relevant seg table
  seg_list <- segobj$segtab
  Chromosome = seg_list[,"Chromosome"] 
  seg_list = as.matrix( seg_list[,-1] )  ## HACK - take out non-numeric Chromosome field to fix below
  
  # make vectors of 0s for columns
  modal_cn = vector(mode="numeric", length=nrow(seg_list))
  expected_cn = vector(mode="numeric", length=nrow(seg_list))
  hz = vector(mode="numeric", length=nrow(seg_list))
  subclonal = vector(mode="numeric", length=nrow(seg_list))
  copy_ratio = vector(mode="numeric", length=nrow(seg_list))

#  cancer_cell_frac = vector(mode="numeric", length=nrow(seg_list))
#  ccf_ci95_low =  vector(mode="numeric", length=nrow(seg_list))
#  ccf_ci95_high =  vector(mode="numeric", length=nrow(seg_list))
  
#  sc_tab = segobj$mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[1, , ]
#  ccf_hat = round(sc_tab[, "CCF_hat"], 5 )
#  ccf_ci95= round(sc_tab[, c("CI95_low", "CI95_high")], 5)

  # print modal and expected values and check for hz and LOH
  for (i in seq_len(nrow(seg_list))) {
    cn = seg_list[i, "copy_num" ]
    copy_ratio[i] = round(sum(cn) / 2, 5)
    modal_cn[i] = max_mat[i] - 1
    expected_cn[i] = round(exp_mat[i] - 1, 5)
#    subclonal[i] = subclonal_ix[i]
#    cancer_cell_frac[i] = ccf_hat[i]
#    ccf_ci95_low[i] = ccf_ci95[i, 1]
#    ccf_ci95_high[i] = ccf_ci95[i, 2]
    
    if (modal_cn[i] == 0) {
      hz[i] = 1
    }
  }
  
  # round and delete appropriate fields from existing seg table
  ix = which(colnames(seg_list) %in% c("copy_num"))

  tab = round(seg_list[, c(-ix)], 5)
## Add back text Chromosome col
  tab = data.frame( Chromosome, tab, stringsAsFactors=FALSE )
    
  return(cbind(tab, copy_ratio, modal_cn, expected_cn, subclonal, cancer_cell_frac, 
               ccf_ci95_low, ccf_ci95_high, hz))
}

