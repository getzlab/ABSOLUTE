## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

AllelicCalcChrArmDistr <- function(seg.obj, seg.q, chr.arms.dat) {
  n.arm <- nrow(chr.arms.dat)
  chr.arm.tab <- array(NA, dim=c(2, n.arm, ncol(seg.q)))
  
  for (i in 1:n.arm) {
    chr.dat <- GetAllelicChrArmSegs(seg.obj, chr.arms.dat[i, ])      
    
    if (length(chr.dat[["int.w"]]) == 0 ) {
      next
    }
    
    chr.arm <- array(NA, dim=c(2, ncol(seg.q)))
    
    chr.arm[1, ] <- colSums(seg.q[chr.dat[["low.ix"]], , drop=FALSE] *
      chr.dat[["int.w"]])
    chr.arm[2, ] <- colSums(seg.q[chr.dat[["high.ix"]], , drop=FALSE] *
      chr.dat[["int.w"]])
    
    chr.arm.tab[, i, ] <- 0
    chr.arm.tab[1, i, which.max(chr.arm[1,])] <- 1
    chr.arm.tab[2, i, which.max(chr.arm[2,])] <- 1
  }  
  
  return(chr.arm.tab)
}


AllelicGetCopyRatioComb <- function(q, delta, b, hscn.params) {
    xx <- (2 * delta * (c(1:q) - 1) + b)
    return(xx)
}

allelic_get_cr_grid_from_ccf_grid = function(qc, qs, delta, comb, ccf_grid) {
  d = qs - qc
  cr_grid = (2 * (delta * d) * ccf_grid) + comb[qc + 1]
  
  if (d < 0) {
    cr_grid = rev(cr_grid)
  }
  
  return(cr_grid)
}

GetAllelicChrArmSegs <- function(seg.obj, chr.arm.dat) {
    seg.dat <- seg.obj[["as.seg.dat"]]
    
    ##  chr, start, stop   
    uniq.seg.ix <- unique(seg.dat[, "seg.ix"])
    seg.ix <- seg.dat[, "seg.ix"]
    n.seg <- length(uniq.seg.ix)
    
    int.w <- c()
    low.ix <- c()
    high.ix <- c()
    
    arm.len.bp <- chr.arm.dat["End.bp"] - chr.arm.dat["Start.bp"]
    
    for (i in 1:n.seg) {
        if (sum(seg.ix == uniq.seg.ix[i]) != 2) {
            next
        }
        
        ix <- which(seg.ix == uniq.seg.ix[i])[1]
        
        if (seg.dat[ix, "Chromosome"] != chr.arm.dat["chr"]) {
            next
        }
        
        int.start <- max(as.numeric(c(seg.dat[ix, "Start.bp"],
                                      chr.arm.dat["Start.bp"])))
        int.end <- min(as.numeric(c(seg.dat[ix, "End.bp"],
                                    chr.arm.dat["End.bp"])))
        
        if (int.start > int.end)  {
          next
        }  ## seg does not overlap region
        
        pair.ix <- which(seg.ix == uniq.seg.ix[i])
        
        int.len.bp <- int.end - int.start
        int.w <- c(int.w, int.len.bp / arm.len.bp)
        
        low.ix <- c(low.ix, pair.ix[which.min(seg.dat[pair.ix, "copy_num"])])
        high.ix <- c(high.ix, pair.ix[which.max(seg.dat[pair.ix, "copy_num"])])
    }
    
    int.w <- as.numeric(int.w)
    
    return(list(int.w=int.w, low.ix=low.ix, high.ix=high.ix))
}

AllelicGetAbsSegDat <- function(segobj) 
{
   as.segs = get_allelic_abs_CN_seg_dat_from_allelic_CAPSEG_obj(segobj) 
   tot.segs = get_total_abs_CN_seg_dat_from_allelic_CAPSEG_obj(segobj)

## Graft two tabs together with shared colunms in front.
   shared.cols = intersect( colnames(tot.segs), colnames(as.segs))

## Merge in total CN segs that are absent in allelic tab
   tot.chrend = data.frame(tot.segs[,c("Chromosome", "End.bp")])
   as.chrend = data.frame( as.segs[,c("Chromosome", "End.bp")] )
 
# Assume as segs were strictly filtered from set of total segs, with no change in any breakpoints.
# Assume that all as segs are also total cn segs
   tot_keys = paste( tot.chrend[,1], tot.chrend[,2], "__" )
   as_keys = paste( as.chrend[,1], as.chrend[,2], "__" )
   m.ix = match( tot_keys, as_keys )

   tot.only.ix = is.na(m.ix)

   allelic.cols = setdiff(colnames(as.segs), shared.cols)
   merged = data.frame( tot.segs,  matrix( NA, nrow=nrow(tot.segs), ncol=length(allelic.cols)))
   colnames(merged) = c( colnames(tot.segs), allelic.cols)

# fill in allelic cols where available
   merged[ !tot.only.ix, allelic.cols] = as.segs[ m.ix[!tot.only.ix], allelic.cols ] 

#  return(cbind(tab, "total_copy_ratio"=copy_ratio, "modal_total_cn"=modal_cn, "expected_total_cn"=expected_cn, "total_HZ"=HZ, "total_amp"=amp, "corrected_total_cn"=rescaled_cn) )

# overwrite total CN columns with allelic sums, where available
#   merged[ !tot.only.ix, "corrected_total_cn" ] = rowSums(merged[ !tot.only.ix, c("rescaled.cn.a1", "rescaled.cn.a2") ])

   return( merged )
} 

## For output of total CN ABS model
get_total_abs_CN_seg_dat_from_allelic_CAPSEG_obj = function(segobj)
{
 # Get column number of the max of each row and the expected 
## ONLY USE THESE FIELDS FROM MODEL: 

  seg.amp <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["tot.seg.amp.tab"]]
  seg_qz_tab <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["tot.seg.qz.tab"]]
  seg_q_tab <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["tot.seg.q.tab"]]
##
  
  Q = ncol(seg_q_tab)
  qq = Q

  max_mat <- apply(seg_qz_tab, MARGIN=1, function(x) which.max(x))
  subclonal_ix = (max_mat == (Q + 1))
  
  max_mat = apply(seg_q_tab, MARGIN=1, function(x) which.max(x))
  
  exp_mat = apply(seg_q_tab, MARGIN=1, 
                  function(x) {
                    x <- x[1:qq] / sum(x[1:qq])
                    return(sum(x * c(1:qq))) 
                  })
  
  # seg_list is relevant seg table
  seg_list <- segobj$total.seg.dat
  Chromosome = seg_list[,"Chromosome"] 
  seg_list = as.matrix( seg_list[,-1] )  ## HACK - take out non-numeric Chromosome field to fix below
  
  # make vectors of 0s for columns
  modal_cn = vector(mode="numeric", length=nrow(seg_list))
  expected_cn = vector(mode="numeric", length=nrow(seg_list))
  HZ = vector(mode="numeric", length=nrow(seg_list))
  subclonal = vector(mode="numeric", length=nrow(seg_list))
  amp <- vector(mode="numeric", length=nrow(seg_list))
  copy_ratio = vector(mode="numeric", length=nrow(seg_list))
  corrected_total_cn = vector(mode="numeric", length=nrow(seg_list))
  rescaled_total_cn = vector(mode="numeric", length=nrow(seg_list))


#  cancer_cell_frac = vector(mode="numeric", length=nrow(seg_list))
#  ccf_ci95_low =  vector(mode="numeric", length=nrow(seg_list))
#  ccf_ci95_high =  vector(mode="numeric", length=nrow(seg_list))
  
#  sc_tab = segobj$mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[1, , ]
#  ccf_hat = round(sc_tab[, "CCF_hat"], 5 )
#  ccf_ci95= round(sc_tab[, c("CI95_low", "CI95_high")], 5)

  alpha = segobj[["mode.res"]][["mode.tab"]][1,"alpha"]
  tau = segobj[["mode.res"]][["mode.tab"]][1,"tau"]

  res = get_b_and_delta( alpha, tau) 
  b= res$b
  delta = res$delta

  # print modal and expected values and check for HZ and LOH
  for (i in seq_len(nrow(seg_list))) 
  {
    cn = seg_list[i, "copy_num" ]
    copy_ratio[i] = round(sum(cn) / 2, 5)
    modal_cn[i] = max_mat[i] - 1
    expected_cn[i] = round(exp_mat[i] - 1, 5)
    amp[i] = seg.amp[i]   


## rescale amps (off comb) by purity / ploidy 
    if( amp[i] )
    {
       corrected_total_cn[i] = (copy_ratio[i] - b) / delta
    }
    else{ corrected_total_cn[i] = expected_cn[i] }

    rescaled_total_cn[i] = (copy_ratio[i] - b) / delta
#    subclonal[i] = subclonal_ix[i]
#    cancer_cell_frac[i] = ccf_hat[i]
#    ccf_ci95_low[i] = ccf_ci95[i, 1]
#    ccf_ci95_high[i] = ccf_ci95[i, 2]
    
    if (modal_cn[i] == 0) {
      HZ[i] = 1
    }
  }
  
  # round and delete appropriate fields from existing seg table
  ix = which(colnames(seg_list) %in% c("copy_num"))
  tab = round(seg_list[, c(-ix)], 5)
## Add back text Chromosome col
  tab = data.frame( Chromosome, tab, stringsAsFactors=FALSE )

## renormalize seg genome fractions
  tab[,"W"] = tab[,"W"] / sum(tab[,"W"])
    
  return(cbind(tab, "total_copy_ratio"=copy_ratio, "modal_total_cn"=modal_cn, "expected_total_cn"=expected_cn, "total_HZ"=HZ, "total_amp"=amp, "corrected_total_cn"=corrected_total_cn, "rescaled_total_cn"=rescaled_total_cn) )
}
   


get_allelic_abs_CN_seg_dat_from_allelic_CAPSEG_obj <- function(segobj) {
  
  ## Get column number of the max of each row and the expected 
  seg.amp <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["seg.amp.tab"]]
#  seg.qz.tab <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["seg.qz.tab"]]
  seg.q.tab <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["seg.q.tab"]]

  Q = ncol(seg.q.tab)
  qq <- Q
  
  max.mat <- apply(seg.q.tab,1, which.max)
  exp.mat <- apply(seg.q.tab, 1, function(x) {
    x <- x[1:Q] / sum(x[1:Q])
    return(sum(x * c(1:Q)))
  })
  
  ## seg_list is relevant seg table
  seg.list <- segobj[["as.seg.dat"]]
  Chromosome = seg.list[,"Chromosome"] 
  seg.list = as.matrix( seg.list[,-1] )  ## HACK - take out non-numeric Chromosome field to fix below
 
  ## get unique seg_ix values
  u <- unique(seg.list[, "seg.ix"])
  
  ## make vectors of 0s for columns
  modal.a1 <- vector(mode="numeric", length=length(u))
  expected.a1 <- vector(mode="numeric", length=length(u))
  modal.a2 <- vector(mode="numeric", length=length(u))
  expected.a2 <- vector(mode="numeric", length=length(u))
  LOH <- vector(mode="numeric", length=length(u))
  HZ <- vector(mode="numeric", length=length(u))
  SC_HZ <- vector(mode="numeric", length=length(u))
  subclonal.a1 <- vector(mode="numeric", length=length(u))
  subclonal.a2 <- vector(mode="numeric", length=length(u))
  sc.qs.a1 <- vector(mode="numeric", length=length(u))
  sc.qs.a2 <- vector(mode="numeric", length=length(u))

  amp.a1 <- vector(mode="numeric", length=length(u))
  amp.a2 <- vector(mode="numeric", length=length(u))

  rescaled.cn.a1 <- vector(mode="numeric", length=length(u))
  rescaled.cn.a2 <- vector(mode="numeric", length=length(u))

  copy.ratio <- vector(mode="numeric", length=length(u))
  hscr.a1 <- vector(mode="numeric", length=length(u))
  hscr.a2 <- vector(mode="numeric", length=length(u))
  
  cancer.cell.frac.a1 <- vector(mode="numeric",length=length(u))
  cancer.cell.frac.a2 <- vector(mode="numeric",length=length(u))
  
  ccf.ci95.low.a1 <-  vector(mode="numeric",length=length(u))
  ccf.ci95.low.a2 <-  vector(mode="numeric",length=length(u))
  
  ccf.ci95.high.a1 <-  vector(mode="numeric",length=length(u))
  ccf.ci95.high.a2 <-  vector(mode="numeric",length=length(u))
  
  sc_tab = segobj[["mode.res"]][["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]][1,,]
  subclonal.ix = sc_tab[,"subclonal_ix"]

  ccf_hat = round(sc_tab[, "CCF_hat"], 5)
  ccf_ci95 = round(sc_tab[, c("CI95_low", "CI95_high")], 5)
  sc_qs = sc_tab[,"qs"]
  
  alpha = segobj[["mode.res"]][["mode.tab"]][1,"alpha"]
  tau = segobj[["mode.res"]][["mode.tab"]][1,"tau"]
  res = get_b_and_delta( alpha, tau) 
  b= res$b
  delta = res$delta

  ## for each pair in seg_ix, print modal and expected values and check for HZ and LOH
  for (i in seq_along(u)) {
    seg <- u[i]
    usegs <- which(seg.list[, "seg.ix"] == seg)
    hscr <- sort(seg.list[usegs, "copy_num"])
    copy.ratio[i] <- round(sum(hscr) / 2, 5)
    hscr.a1[i] <- round(hscr[1], 5)
    hscr.a2[i] <- round(hscr[2], 5)
    
    modal.a1[i] <- max.mat[usegs[1]] - 1
    modal.a2[i] <- max.mat[usegs[2]] - 1
    expected.a1[i] <- round(exp.mat[usegs[1]] - 1, 5)
    expected.a2[i] <- round(exp.mat[usegs[2]] - 1, 5)
    
    subclonal.a1[i] <- subclonal.ix[usegs[1]]
    subclonal.a2[i] <- subclonal.ix[usegs[2]]
    sc.qs.a1[i] = sc_qs[usegs[1]]
    sc.qs.a2[i] = sc_qs[usegs[2]]
    
    amp.a1[i] = seg.amp[usegs[1]]   
    amp.a2[i] = seg.amp[usegs[2]]   


## rescale amps (off comb) by purity / ploidy 
    if( amp.a1[i] )
    {
       rescaled.cn.a1[i] = (hscr.a1[i] - b) / (2*delta)
    }
    else{ rescaled.cn.a1[i] = expected.a1[i] }

    if( amp.a2[i] )
    {
       rescaled.cn.a2[i] = (hscr.a2[i] - b) / (2*delta)
    }
    else{ rescaled.cn.a2[i] = expected.a2[i] }


    cancer.cell.frac.a1[i] = ccf_hat[usegs[1]]
    cancer.cell.frac.a2[i] = ccf_hat[usegs[2]]
    ccf.ci95.low.a1[i] = ccf_ci95[usegs[1], 1]
    ccf.ci95.high.a1[i] = ccf_ci95[usegs[1], 2]
    ccf.ci95.low.a2[i] = ccf_ci95[usegs[2], 1]
    ccf.ci95.high.a2[i] = ccf_ci95[usegs[2], 2]
    
    if ((modal.a1[i]==0) && (modal.a2[i]==0)) {
      HZ[i] <- 1
    }

## This is not safe since sc.qs.a. can be NA when neg.ix, i.e., when way below 0
#    else if ( (modal.a2[i] == 0 & subclonal.a1[i]==1 & sc.qs.a1[i]==0) |
#              (modal.a1[i] == 0 & subclonal.a2[i]==1 & sc.qs.a2[i]==0) )
#    {
#      SC_HZ[i] = 1
#    } 

    if ((modal.a1[i] == 0) || (modal.a2[i]==0)) {
      LOH[i] <- 1
    }
  }       
  
  ## round and delete appropriate fields from existing seg table
  ix <- which(colnames(seg.list) %in% c("seg.ix", "copy_num"))
  tab <- round(seg.list[1:length(u), c(-ix)], 5)



## Add back text Chromosome col
  u_chr = Chromosome[1:length(u)]
  tab = data.frame( "Chromosome"=u_chr, tab, stringsAsFactors=FALSE, check.names=FALSE )
   
  res = data.frame( tab, copy.ratio, hscr.a1, hscr.a2, modal.a1, modal.a2,
               expected.a1, expected.a2, subclonal.a1, subclonal.a2,
               cancer.cell.frac.a1, ccf.ci95.low.a1, ccf.ci95.high.a1,
               cancer.cell.frac.a2, ccf.ci95.low.a2, ccf.ci95.high.a2,
               LOH, HZ, SC_HZ, amp.a1, amp.a2, rescaled.cn.a1, rescaled.cn.a2,
		stringsAsFactors=FALSE, check.names=FALSE )

  return(res)
}
