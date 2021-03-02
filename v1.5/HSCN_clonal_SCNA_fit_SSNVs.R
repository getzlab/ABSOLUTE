## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GetHscnSomaticMutComb <- function(alpha, q1, q2) {
  q.a <- c(1:q2)
  q <- q1 + q2
  f.reads <- alpha * q.a/(2 * (1 - alpha) + alpha * q)
  vals <- q.a
  res <- cbind(vals, f.reads)
  return(res)
}


GetHscnSomaticXMutComb <- function(alpha, q1, q2) {
  q.a <- c(1:q2)
  q <- q1 + q2
  f.reads <- alpha * q.a / (1 * (1 - alpha) + alpha * q)
  vals <- q.a
  res <- cbind(vals, f.reads)
  return(res)
}

# eqn. 11
GetHscnGermlineMutComb <- function(alpha, q1, q2) {
  eps <- 0.001
  
  q.vals <- c(q1, q2, q1 + q2)
  Q <- q1 + q2
  f.reads <- ((1 - alpha) + alpha * q.vals)/(2 * (1 - alpha) + alpha * Q)
  f.reads[3] <- 1 - eps
  vals <- q.vals
  res <- cbind(vals, f.reads)
  
  return(res)
}

GetHscnGermlineXMutComb <- function(alpha, q1, q2) {
  eps <- 0.001
  
  q.vals <- c(q1, q2, q1 + q2)
  Q <- q1 + q2
  f.reads <- ((1 - alpha) + alpha * q.vals)/((1 - alpha) + alpha * Q)
  f.reads[3] <- 1 - eps
  vals <- q.vals
  res <- cbind(vals, f.reads)
  
  return(res)
}



AllelicGetMutSegIx <- function(maf, segtab) {
  ## compute lookup table for each mutation into seg.q.tab
  mut.seg.ix <- matrix(NA, nrow=nrow(maf), ncol=2)
  colnames(mut.seg.ix) <- c("A1.ix", "A2.ix")

  ## This is Dan-Avi's fault
  if (!("End_position" %in% colnames(maf))) {
    end.position <- maf[, "Start_position"]
    maf <- cbind(maf, "End_position"=end.position)
  }
  
  for (i in seq_len(nrow(maf)))
  {
     seg.ix <- (maf[ i, "Chromosome" ] == segtab[,"Chromosome"]) &
               (maf[ i, "Start_position"] >= segtab[,"Start.bp"]) &
               (maf[ i, "End_position"] <= segtab[,"End.bp"])

## not found on both homologues?
     if (sum(seg.ix) != 2) {
       next
     }
    
     seg.ids <- which(seg.ix)
     mut.seg.ix[i, 1] <- seg.ids[1] 
     mut.seg.ix[i, 2] <- seg.ids[2] 
  }

## try to map missing muts
   nix <- which(apply(is.na(mut.seg.ix), 1, sum) > 0)
   if( length(nix) > 0 )
   {
      for( i in 1:length(nix))
      {
         m.ix = nix[i]
         if( !(maf[m.ix,"Chromosome"] %in% segtab[,"Chromosome"])) { next }  #no contig

   ## is mut closer to end or current seg?, or start of next seg?
   ## ASSUME genomic ordering of segs!!

     # find s.ix = seg before breakpoint with SSNV inside
         chr.ix = segtab[,"Chromosome"] == maf[m.ix,"Chromosome"]
         chrtab = segtab[chr.ix,]

         s.ix = which.min( abs(maf[m.ix,"Start_position"] - chrtab[,"End.bp"]) )

         s.len = maf[m.ix, "Start_position"] - chrtab[s.ix,"End.bp"] 
         e.len = chrtab[s.ix+1,"Start.bp"] - maf[m.ix, "End_position"]

         if( s.len <= e.len ) { use.ix = s.ix }
         if( s.len > e.len ) { use.ix = s.ix+1 }

         seg.ix = 	chrtab[ use.ix, "Chromosome"] == segtab[,"Chromosome"] & 
			chrtab[ use.ix, "Start.bp"] == segtab[,"Start.bp"] & 
			chrtab[ use.ix, "End.bp"] == segtab[,"End.bp"]  

## not found on both homologues?
         if (sum(seg.ix) != 2) {
           next
         }
    
         seg.ids <- which(seg.ix)
         mut.seg.ix[m.ix, 1] <- seg.ids[1] 
         mut.seg.ix[m.ix, 2] <- seg.ids[2] 
      }
   }

  return(mut.seg.ix)
}




get_hap.q.keys = function(mut.cn.dat)
{
  hap.q.keys <- paste(mut.cn.dat[, "HS_q_hat_1"], mut.cn.dat[, "HS_q_hat_2"],
                      sep="__" )
  
  u.keys<- unique(hap.q.keys)
  n.keys <- length(u.keys)
  mut.n.keys <- rep(NA, n.keys)
  
  for(i in seq_len(n.keys)) {
    res <- strsplit(u.keys[i], "__" )[[1]]
    q1 <- as.integer(res[1])
    q2 <- as.integer(res[2])
    mut.ix <- (mut.cn.dat[, "HS_q_hat_1"] == q1) & (mut.cn.dat[, "HS_q_hat_2"] == q2)

    mut.n.keys[i] <- sum(mut.ix)
  }
  
  u.keys <- u.keys[ order(mut.n.keys, decreasing=TRUE) ]
  mut.n.keys <- sort(mut.n.keys, decreasing=TRUE)

  return(list(u.keys=u.keys, mut.n.keys=mut.n.keys, n.keys=n.keys))
}




AllelicCalcGermlineSNVLoglik <- function(mut.cn.dat, mode_info, SSNV_model)
{
  mut_class_w = SSNV_model[["mut_class_w"]]

  alpha = mode_info["alpha"]
  keys = get_hap.q.keys(mut.cn.dat)
  u.keys=keys$u.keys; mut.n.keys=keys$mut.n.keys; n.keys=keys$n.keys

  LL_mat = matrix(NA, nrow=nrow(mut.cn.dat), ncol=3 )
 
  for (i in seq_len(n.keys)) 
  {
    res <- strsplit(u.keys[i], "__" )[[1]]
    q1 <- as.integer(res[1])
    q2 <- as.integer(res[2])
    mut.ix <- (mut.cn.dat[, "HS_q_hat_1"] == q1) & (mut.cn.dat[, "HS_q_hat_2"] == q2) 

    ## Germline model calc
    log.gl.w <- log(GermlineMutPrior(mut_class_w))
    
    germ.comb <- GetHscnGermlineMutComb(alpha, q1, q2)
    log.gl.w <-  matrix(log.gl.w, ncol=nrow(germ.comb), nrow=length(alt), byrow=TRUE) 
    LL_mat[mut.ix,] <- log.gl.w + SSNV_FhatCombPost(alt, ref, germ.comb[,"f.reads"], SSNV_model )
  }

  return(LL_mat)
}




AllelicCalcClonalSSNVLoglik <- function(mut.cn.dat, mode_info, SSNV_model)
{
  som.theta.q = SSNV_model[["som_theta_Q_mode"]]
  mut_class_w = SSNV_model[["mut_class_w"]]

  LL_mat <- matrix(-Inf, nrow=nrow(mut.cn.dat), ncol=length(som.theta.q))

  alpha = mode_info["alpha"]
  keys = get_hap.q.keys(mut.cn.dat)
  u.keys=keys$u.keys; mut.n.keys=keys$mut.n.keys; n.keys=keys$n.keys
 
  for (i in seq_len(n.keys)) 
  {
    res <- strsplit(u.keys[i], "__" )[[1]]
    q1 <- as.integer(res[1])
    q2 <- as.integer(res[2])
    mut.ix <- (mut.cn.dat[, "HS_q_hat_1"] == q1) & (mut.cn.dat[, "HS_q_hat_2"] == q2) 

    ## FIXME: Is this really a &&?
    if ((q1 == 0) & (q2 == 0)) {
      LL_mat[mut.ix, ] <- -Inf
      next
    }

    alt = mut.cn.dat[mut.ix, "alt"]
    ref =  mut.cn.dat[mut.ix, "ref"]

## TODO: add support for X-chr
#    m.ix = mut.cn.dat[mut.ix, "male_X"]

    ## Somatic autosome model calc
    som.comb <- GetHscnSomaticMutComb(alpha, q1, q2)
    som.w <- clonal_SSNV_Mult_Prior(mut_class_w, nrow(som.comb), som.theta.q)
    som.w <-  matrix(som.w, ncol=nrow(som.comb), nrow=length(alt), byrow=TRUE)
    som.ll <- log(som.w) + SSNV_FhatCombPost( alt, ref, som.comb[,"f.reads"], SSNV_model )

    LL_mat[mut.ix, c(1:ncol(som.ll))] = som.ll
  }

  return(LL_mat)
}   


## Replaces: "CalcSampleMutsPostPr"
allelic_eval_SNV_models_evidence = function(mut.cn.dat, mode_info, post.ccf.LL.grid, SSNV_model )
{
  ## n_mut X Q 
   clonal.comb.LL = AllelicCalcClonalSSNVLoglik(mut.cn.dat, mode_info, SSNV_model)
   Z = LogAdd(clonal.comb.LL)
   Z[is.nan(Z)] = -Inf
   som_mut_Q_tab = exp(clonal.comb.LL - Z)
   homdel.ix = mut.cn.dat[, "HS_q_hat_1"]==0 & mut.cn.dat[, "HS_q_hat_2"] == 0  

  ## n_mut X 1
   som.clonal.log.ev = LogAdd(clonal.comb.LL) #+ log(SSNV_model[["mut_class_w"]][["SM"]])
   som.clonal.log.ev[is.nan(som.clonal.log.ev)] = -Inf

#   subclonal.exp.log.ev = H_exp_subclonal_SSNV( post.ccf.LL.grid, SSNV_model )
   subclonal.exp.log.ev = rep(-Inf, nrow(mut.cn.dat))
   subclonal.unif.log.ev = H_unif_subclonal_SSNV( post.ccf.LL.grid, SSNV_model )

   if( SSNV_model[["mut_class_w"]][["GL"]] > 0 )
   {
      germline.comb.LL = AllelicCalcGermlineSNVLoglik(mut.cn.dat, mode_info, SSNV_model)
      germline.log.ev = LogAdd(germline.comb.LL) #+ log(SSNV_model[["mut_class_w"]][["GL"]])
      colnames(germline.log.ev) = c("alt_minor", "alt_major", "hom_alt")
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

   post_Prs = allelic_annotate_SSNVs_on_clonal_SCNAs( mut.ev.mat, som_mut_Q_tab, mut.cn.dat, SSNV_model, LL)

   return(list("mut.ev.mat"=mut.ev.mat, "som_mut_Q_tab"=som_mut_Q_tab, "post_Prs"=post_Prs))
}


allelic_annotate_SSNVs_on_clonal_SCNAs = function( mut.ev.mat, som_mut_Q_tab, mut.cn.dat, SSNV_model, LL )
{
   mut.w = SSNV_model[["mut_class_w"]]

# Original contract for CalcSampleMutsPostPr()
   cols <- c("Pr_somatic_clonal", "Pr_germline", "Pr_subclonal",
             "Pr_subclonal_wt0", "Pr_wt0", "Pr_ge2", "Pr_GL_som_HZ_alt",
             "Pr_GL_som_HZ_ref", "Pr_cryptic_SCNA", "modal_q_s", "LL")
   post_Prs = matrix( 0, nrow=nrow(mut.cn.dat), ncol=length(cols) )
   colnames(post_Prs)=cols
   post_Prs[,"LL"] = LL

   LOH.ix = mut.cn.dat[, "HS_q_hat_1"]==0 & mut.cn.dat[, "HS_q_hat_2"] > 0 
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
      q2 = mut.cn.dat[LOH.ix, "HS_q_hat_2"]
      post_Prs[LOH.ix,"Pr_wt0"] = t(s1)[ q2 + ncol(s1) * r ]  * post_Prs[LOH.ix, "Pr_somatic_clonal"] 
   }

   ix.01 = (mut.cn.dat[, "HS_q_hat_1"]==0 & mut.cn.dat[, "HS_q_hat_2"] == 1)  
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



