## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

allelic_get_muts_nearest_clonal_scna <- function(mut.cn.dat, seg.q.tab, Q) {
  if (!all(c("A1.ix", "A2.ix") %in% colnames(mut.cn.dat)) ) {
    stop("wrong colnames")
  }
  
  muts.p.q <- array(NA, dim=c( nrow(mut.cn.dat), 2, Q))
  
  for (i in seq_len(nrow(mut.cn.dat))) {
    muts.p.q[i, 1, ] <- seg.q.tab[mut.cn.dat[i,"A1.ix"], ]   
    muts.p.q[i, 2, ] <- seg.q.tab[mut.cn.dat[i,"A2.ix"], ]   
  }
  
  ## todo: integrate over this instead
  muts.q.hat <- cbind(apply(muts.p.q[, 1, ], 1, which.max ) - 1,
                      apply(muts.p.q[, 2, ], 1, which.max ) - 1)

  mut.cn.tab <- cbind("q_hat" = rowSums(muts.q.hat),
                      "HS_q_hat_1"=muts.q.hat[,1],
                      "HS_q_hat_2"=muts.q.hat[,2] ) 
  
  return(mut.cn.tab) 
}


allelic_get_subclonal_scna_mut_ix = function(mut.cn.dat, subclonal_scna_tab) 
{
  wix1 = which(subclonal_scna_tab[, "subclonal_ix"] == 1)
  wix2 = which( !is.na(subclonal_scna_tab[, "qs"]) )

  ix = (mut.cn.dat[, "A1.ix"] %in% wix1 | mut.cn.dat[, "A2.ix"] %in% wix1 ) &
       (mut.cn.dat[, "A1.ix"] %in% wix2 & mut.cn.dat[, "A2.ix"] %in% wix2 ) 
  
  return(ix)
}


get_allelic_subclonal_scna_tab = function(mut.cn.dat, subclonal_scna_tab, seg.q.tab ) 
{
  # loop over HSCR segs
  a_ix_keys = paste(mut.cn.dat[, "A1.ix"], mut.cn.dat[, "A2.ix"], sep="__" )
  u_keys= unique(a_ix_keys)
  n_keys = length(u_keys)
  mut_n_keys = rep(NA, n_keys)

  tabcols = c("SC.ix", "C.ix", "total_qc", "total_qs", "SC_Aq_d", "SC_Aq_a", "C_Aq", "A1.ix", "A2.ix")
  allelic_subclonal_scna_tab = matrix( NA, nrow=length(u_keys), ncol=length(tabcols) )
  colnames(allelic_subclonal_scna_tab) = tabcols

  if( n_keys == 0 )  ## No muts on subclonal SCNAs!
  {
    subclonal.mut.tab = as.data.frame(matrix(NA, nrow=0, ncol=ncol(allelic_subclonal_scna_tab) ) )
   colnames(subclonal.mut.tab) = colnames(allelic_subclonal_scna_tab)
   subclonal.mut.tab = subclonal.mut.tab[ , setdiff(colnames(subclonal.mut.tab), c("A1.ix","A2.ix"))]
   return( list("subclonal_scna_tab" = allelic_subclonal_scna_tab, "subclonal.mut.tab"=subclonal.mut.tab) )
  }
  
  for (i in 1:n_keys) {
    res = strsplit(u_keys[i], "__" )[[1]]
    A1.ix = as.integer(res[1])
    A2.ix = as.integer(res[2])
    mut_ix = mut.cn.dat[, "A1.ix"] == A1.ix & mut.cn.dat[,"A2.ix"] == A2.ix 
    mut_n_keys[i] = sum(mut_ix)
  }
  
  u_keys = u_keys[order(mut_n_keys, decreasing=TRUE)]
  mut_n_keys = sort(mut_n_keys, decreasing=TRUE)
  rownames(allelic_subclonal_scna_tab) = u_keys

  SC_Prs = matrix(NA, nrow=length(u_keys), ncol=2)
  colnames(SC_Prs) = c("Pr_SC_SC.ix", "Pr_SC_C.ix")

  for( i in 1:n_keys )
  {
    res = strsplit(u_keys[i], "__" )[[1]]
    A1.ix = as.integer(res[1])
    A2.ix = as.integer(res[2])

    allelic_subclonal_scna_tab[i,"A1.ix"] = A1.ix
    allelic_subclonal_scna_tab[i,"A2.ix"] = A2.ix

  ## Choose which homologous seg is subclonal 
# subclonal_scna_tab = cbind(ccf_sum, subclonal_ix, Pr_subclonal = seg.post.subclonal,  
#				qs=CN_states[,"qs"], qc=CN_states[,"qc"] )

    if( subclonal_scna_tab[ A1.ix, "Pr_subclonal"] < subclonal_scna_tab[ A2.ix, "Pr_subclonal"] )
    {
       SC.ix = A2.ix
        C.ix = A1.ix
    }
    else
    {
       SC.ix = A1.ix
        C.ix = A2.ix 
    }

    clonal_allelic_CN = which.max(seg.q.tab[C.ix,]) - 1

    allelic_subclonal_scna_tab[i,"SC.ix"] = SC.ix
    allelic_subclonal_scna_tab[i,"C.ix"] = C.ix

    SC_Prs[i,"Pr_SC_SC.ix"] = subclonal_scna_tab[ SC.ix, "Pr_subclonal" ]
    SC_Prs[i,"Pr_SC_C.ix"] = subclonal_scna_tab[ C.ix, "Pr_subclonal" ]

    allelic_subclonal_scna_tab[i, "total_qc"] = clonal_allelic_CN + subclonal_scna_tab[ SC.ix, "qc"]
    allelic_subclonal_scna_tab[i, "total_qs"] = clonal_allelic_CN + subclonal_scna_tab[ SC.ix, "qs"]
    allelic_subclonal_scna_tab[i, "SC_Aq_d"] = subclonal_scna_tab[ SC.ix, "qs"]
    allelic_subclonal_scna_tab[i, "SC_Aq_a"] = subclonal_scna_tab[ SC.ix, "qc" ]
    allelic_subclonal_scna_tab[i, "C_Aq"] = clonal_allelic_CN
  }

  allelic_subclonal_scna_tab = data.frame( allelic_subclonal_scna_tab, SC_Prs, check.names=FALSE, stringsAsFactors=FALSE )

  if( any(is.na(allelic_subclonal_scna_tab))) { stop("NA in allelic_subclonal_scna_tab") }

  subclonal.mut.tab = as.data.frame(matrix(NA, nrow=nrow(mut.cn.dat), ncol=ncol(allelic_subclonal_scna_tab) ) )
  for( i in 1:nrow(allelic_subclonal_scna_tab) )
  {
    A1.ix = allelic_subclonal_scna_tab[i,"A1.ix"]
    A2.ix = allelic_subclonal_scna_tab[i,"A2.ix"]
    mut_ix = which(mut.cn.dat[,"A1.ix"] == A1.ix & mut.cn.dat[,"A2.ix"] == A2.ix)
    for( j in 1:length(mut_ix))
    {
       subclonal.mut.tab[mut_ix[j],] =  allelic_subclonal_scna_tab[i,]
    }
  }
  colnames(subclonal.mut.tab) = colnames(allelic_subclonal_scna_tab)
  subclonal.mut.tab = subclonal.mut.tab[ , setdiff(colnames(subclonal.mut.tab), c("A1.ix","A2.ix"))]

  if( any(is.na(subclonal.mut.tab))) { stop() }

  return( list("subclonal_scna_tab" = allelic_subclonal_scna_tab, "subclonal.mut.tab"=subclonal.mut.tab) )
}



allelic_calc_sample_muts_on_subclonal_scna = function(mut.cn.dat, mode_info, allelic_subclonal_scna_tab, scna_log_ccf_dens, SSNV_model ) 
{
  kQ = SSNV_model[["kQ"]]
  alpha = mode_info["alpha"]

  N_H = 4
## LL for each mut x H*, assuming mut is clonal (on subclonal SCNA)
  H.123_qm_ll = array( -Inf, dim=c( N_H, nrow(mut.cn.dat), kQ) )

  for( i in 1:nrow(allelic_subclonal_scna_tab) )
  {
    SC.ix = allelic_subclonal_scna_tab[i,"SC.ix"]
    A1.ix = allelic_subclonal_scna_tab[i,"A1.ix"]
    A2.ix = allelic_subclonal_scna_tab[i,"A2.ix"]

    log_f_c_dens = scna_log_ccf_dens[SC.ix, , drop=FALSE]
    mut_ix = mut.cn.dat[,"A1.ix"] == A1.ix & mut.cn.dat[,"A2.ix"] == A2.ix 
    alt = mut.cn.dat[mut_ix, "alt"]
    ref = mut.cn.dat[mut_ix, "ref"]

    res = seg_SSNV_on_subclonal_SCNA_log_ev( alpha, alt, ref, log_f_c_dens, SSNV_model, allelic_subclonal_scna_tab[i,] )

    for( j in 1:(N_H-1) ) {   ## H4 does not have a clonal possibility (LL= -Inf)
       H.123_qm_ll[ j, mut_ix, c(1:ncol(res[[j]]))] = res[[j]]
    }
  }

 ## integrate over qm
  H.123.clonal.log.ev = matrix( apply( H.123_qm_ll, 1, LogAdd ), ncol=N_H, byrow=FALSE )
  H.123.clonal.log.ev[is.nan(H.123.clonal.log.ev)] = -Inf

  H.123_Pr_q = array(NA, dim=dim(H.123_qm_ll) )
  for( i in 1:N_H )
  {
     H.123_Pr_q[i,,] = exp(H.123_qm_ll[i,,] - H.123.clonal.log.ev[,i] )
  }
  H.123_Pr_q[is.nan(H.123_Pr_q)] = 0




  res = calc_SSNV_on_subclonal_SCNA_CCF_dens( alpha, mut.cn.dat, allelic_subclonal_scna_tab, scna_log_ccf_dens, SSNV_model )
  H.123.ssnv.ccf.dens = res[["H.123.ssnv.ccf.dens"]] ## conditional on H* hypotheses
  H.123.ssnv.log.Z = res[["H.123.ssnv.log.Z"]] # log partition func of each H*

  H1.log.ev = LogAdd( cbind( H.123_qm_ll[1,,,drop=FALSE], H.123.ssnv.log.Z[,1,drop=FALSE]) )
  H2.log.ev = LogAdd( cbind( H.123_qm_ll[2,,,drop=FALSE], H.123.ssnv.log.Z[,2,drop=FALSE]) )
  H3.log.ev = LogAdd( cbind( H.123_qm_ll[3,,,drop=FALSE], H.123.ssnv.log.Z[,3,drop=FALSE]) )
  H4.log.ev = H.123.ssnv.log.Z[,4,drop=FALSE]


  H.log.ev = cbind(H1.log.ev, H2.log.ev, H3.log.ev, H4.log.ev) + log(1/4)
  colnames(H.log.ev) = c("H1", "H2", "H3", "H4")

  H.log.ev[is.nan(H.log.ev)] = -Inf
  H.log.Z = LogAdd(H.log.ev)
  H.ev = exp(H.log.ev - H.log.Z)

  ssnv.ccf.dens = H.ev[,1] * H.123.ssnv.ccf.dens[1,,] + 
                  H.ev[,2] * H.123.ssnv.ccf.dens[2,,] +
                  H.ev[,3] * H.123.ssnv.ccf.dens[3,,] +
                  H.ev[,4] * H.123.ssnv.ccf.dens[4,,] 

## TODO: useless - test removing
  subclonal.exp.log.ev = -Inf

## H3 "clonal" model actually implies subclonal SSNV with CCF == SCNA CCF
  subclonal.unif.log.ev = LogAdd( cbind((H.123.ssnv.log.Z + log(H.ev)), H.123.clonal.log.ev[,3] + log(H.ev[,3]) )) + log(SSNV_model[["mut_class_w"]][["SC"]])

  H12.ev = H.ev[,c(1,2)]
  H12.som.clonal.log.ev = H.123.clonal.log.ev[,c(1,2)] + log(H12.ev) + log(SSNV_model[["mut_class_w"]][["SM"]])
  som.clonal.log.ev = LogAdd(H12.som.clonal.log.ev)
  som.clonal.log.ev[is.nan(som.clonal.log.ev)] = -Inf

  if( SSNV_model[["mut_class_w"]][["GL"]] > 0 )
  {
##  TODO: adapt for subclonal SCNA
#     germline.comb.LL = AllelicCalcGermlineSNVLoglik(mut.cn.dat, mode_info, SSNV_model)
#     germline.log.ev = LogAdd(germline.comb.LL) + log(SSNV_model[["mut_class_w"]][["GL"]])
#     colnames(germline.log.ev) = c("alt_minor", "alt_major", "hom_alt")
     stop("Not implemented yet")
  } else {
#     germline.log.ev= matrix(-Inf,ncol=3, nrow=nrow(mut.cn.dat) )
     germline.log.ev= -Inf
  }

  homdel.ix = mut.cn.dat[, "HS_q_hat_1"]==0 & mut.cn.dat[, "HS_q_hat_2"] == 0  

 ## this can't happen in the subclonal SCNA case
  SSNV.on.CN0.log.ev = rep(-Inf, nrow(mut.cn.dat))
  SSNV.on.CN0.log.ev[ homdel.ix ] = SSNV_model[["SSNV.on.CN0"]]

## put together mutually exclusive hypothoses:
  mut.log.ev.mat = cbind( som.clonal.log.ev, subclonal.exp.log.ev, subclonal.unif.log.ev, germline.log.ev, SSNV.on.CN0.log.ev )
  mut.log.ev.mat[homdel.ix, c(1,2,3)] = -Inf 

  LL = LogAdd(mut.log.ev.mat)
  mut.ev.mat = exp(mut.log.ev.mat - LL)

  H12_Pr_clonal = matrix( exp(H12.som.clonal.log.ev - LogAdd( cbind( H12.som.clonal.log.ev, subclonal.unif.log.ev) ) ),  ncol=2  )
## This is over ancestral SSNV multiplicities (H1 and H2 only)
  som_mut_Q_tab = matrix( H12_Pr_clonal[,1] * H.123_Pr_q[1,,] + 
                          H12_Pr_clonal[,2] * H.123_Pr_q[2,,],  
 				ncol=kQ,  nrow=nrow(mut.cn.dat), byrow=FALSE)
## rows sum to rowSums(H12_Pr_clonal)
#  if( any( abs(rowSums(som_mut_Q_tab) - 1) > 1e-10 ) ) { stop() }


  post_Prs = annotate_SSNVs_on_subclonal_SCNAs( mut.ev.mat, H.ev, H.123_Pr_q, som_mut_Q_tab, mut.cn.dat, SSNV_model, H.log.Z)

  result = list( "mut.ev.mat"=mut.ev.mat, "ssnv.ccf.dens"=ssnv.ccf.dens, "H.ev"=H.ev, "som_mut_Q_tab"=som_mut_Q_tab, "post_Prs"=post_Prs, "H.123_Pr_q_a"=H.123_Pr_q, "H.123.ssnv.ccf.dens"=H.123.ssnv.ccf.dens, "H.123.ssnv.log.Z"=H.123.ssnv.log.Z )

  if( any( lapply( lapply(result, is.nan), any ) ) ) { stop("NaN in result") }

  return( result )
}


seg_SSNV_on_subclonal_SCNA_log_ev = function( alpha, alt, ref, log_f_c_dens, SSNV_model, allelic_subclonal_scna_dat )
{
# fc is a vector 
  delta_s = function(fc, qc, qs, alpha) {
    return(alpha / ( 2 * (1 - alpha) + alpha * qc * (1 - fc) + alpha * fc * qs))
  }

  norm_mut_q_w = function(h_qm_ll, som.theta.q) 
  {
 ## Apply and multiplicity (Q) weight
    n = ncol(h_qm_ll)
    qm_w = som.theta.q[seq_len(n)]
      
   # for all muts on seg
    som_w = matrix(qm_w, ncol=n, nrow=nrow(h_qm_ll), byrow=TRUE)
  }

  kQ = SSNV_model[["kQ"]]
  som.theta.q = SSNV_model[["som_theta_Q_mode"]]
  f_c_grid = SSNV_model[["ccf_grid"]]

  tqc = as.integer(allelic_subclonal_scna_dat["total_qc"])
  tqs = as.integer(allelic_subclonal_scna_dat["total_qs"])
  C_Aq = as.integer(allelic_subclonal_scna_dat["C_Aq"])
  SC_Aq_a = as.integer(allelic_subclonal_scna_dat["SC_Aq_a"])
  SC_Aq_d = as.integer(allelic_subclonal_scna_dat["SC_Aq_d"])

## H3: Both SSNV and SCNA are in the same derived cell population:
   H3_max_qm = max(C_Aq, SC_Aq_d)
   if( H3_max_qm > 0 )
   {
     qm_set = seq_len(H3_max_qm)
     int_mat = matrix(log_f_c_dens, nrow=length(qm_set), ncol=length(log_f_c_dens), byrow=TRUE)
     outer_prod = qm_set %*% t( f_c_grid * delta_s(f_c_grid, tqc, tqs, alpha))
     outer_prod[outer_prod>1] = 1  ## rounding error
     h3_qm_ll = mut_qm_fc_grid_integral( alt, ref, int_mat, outer_prod, SSNV_model )
   } else {
     h3_qm_ll = matrix(-Inf, nrow=length(alt), ncol=1)
   }
    
## H2: Ancestral SSNV with derived SCNA in cis with SSNV: (qm can be 0)
   H2_max_qd = SC_Aq_d
   if( SC_Aq_d < SC_Aq_a ) { H2_min_qd = 0 }  # SC loss
   if( SC_Aq_d > SC_Aq_a ) { H2_min_qd = SC_Aq_a }  # SC gain

   qd_set = seq(H2_min_qd, H2_max_qd)  ## derived SSNV mult
   qa_set = qd_set + (tqc-tqs)	       ## ancestral "  "
 ## ancestral multiplicity cannot be 0.
   qd_set = qd_set[qa_set > 0]
   qa_set = qa_set[qa_set > 0]

   if( length(qd_set) > 0 )
   {
      h2_qm_ll = matrix(-Inf, nrow=length(alt), ncol= length(seq(1, max(qa_set))) )
      int_mat = matrix(log_f_c_dens, nrow=length(qd_set), ncol=length(log_f_c_dens), byrow=TRUE)
      outer_prod = (qa_set %*% t((1 - f_c_grid) * delta_s(f_c_grid, tqc, tqs, alpha))) + 
                   (qd_set %*% t(f_c_grid * delta_s(f_c_grid, tqc, tqs, alpha)))  
      outer_prod[outer_prod>1] = 1  ## rounding error

    ## complete deletions in cis cannot be fully clonal 
      outer_prod[ qd_set==0, ncol(outer_prod) ] = 0
      int_mat[ qd_set==0, ncol(outer_prod) ] = -Inf

      h2_qm_ll[,qa_set] = mut_qm_fc_grid_integral( alt, ref, int_mat, outer_prod, SSNV_model )
   } else {
     h2_qm_ll = matrix(-Inf, nrow=length(alt), ncol=1)
   }

## H1: Ancestral SSNV with derived SCNA in trans with SSNV:
   if( SC_Aq_d < SC_Aq_a ) { H1_max_qm = C_Aq }  # SC loss
   if( SC_Aq_d > SC_Aq_a ) { H1_max_qm = SC_Aq_a }  # SC gain
   if(H1_max_qm > 0) {
     qm_set = seq(1, H1_max_qm)
     int_mat = matrix(log_f_c_dens, nrow=length(qm_set), ncol=length(log_f_c_dens), byrow=TRUE)
     outer_prod = qm_set %*% t( delta_s(f_c_grid, tqc, tqs, alpha)) 
     outer_prod[outer_prod>1] = 1  ## rounding error
     h1_qm_ll = mut_qm_fc_grid_integral( alt, ref, int_mat, outer_prod, SSNV_model )
   } else { 
     h1_qm_ll = matrix(-Inf, nrow=length(alt), ncol=1)
   }
 
   som_w = norm_mut_q_w(h1_qm_ll, som.theta.q)
   h1_qm_ll = log(som_w) + h1_qm_ll
   som_w = norm_mut_q_w(h2_qm_ll, som.theta.q)
   h2_qm_ll = log(som_w) + h2_qm_ll 
   som_w = norm_mut_q_w(h3_qm_ll, som.theta.q)
   h3_qm_ll = log(som_w) + h3_qm_ll

#  if( all( !is.finite(cbind(h1_qm_ll, h2_qm_ll, h3_qm_ll)))) { stop() }
   if( any(is.nan(h2_qm_ll))) { stop() }
    
   return(list( "h1_qm_ll"=h1_qm_ll, "h2_qm_ll"=h2_qm_ll, "h3_qm_ll"=h3_qm_ll ))       
}


annotate_SSNVs_on_subclonal_SCNAs = function( mut.ev.mat, H.ev, H.123_Pr_q, som_mut_Q_tab, mut.cn.dat, SSNV_model, LL )
{
   mut.w = SSNV_model[["mut_class_w"]]

# Original contract for CalcSampleMutsPostPr()
   cols <- c("Pr_somatic_clonal", "Pr_germline", "Pr_subclonal",
             "Pr_subclonal_wt0", "Pr_wt0", "Pr_ge2", "Pr_GL_som_HZ_alt",
             "Pr_GL_som_HZ_ref", "Pr_cryptic_SCNA", "modal_q_s", "LL")
   post_Prs = matrix( 0, nrow=nrow(mut.cn.dat), ncol=length(cols) )
   colnames(post_Prs)=cols
   post_Prs[,"LL"] = LL

# mut.cn.dat cols informing SC status:
#  tabcols = c("SC.ix", "C.ix", "total_qc", "total_qs", "SC_Aq_d", "SC_Aq_a", "C_Aq" )
   
   ancestral.LOH.ix = mut.cn.dat[, "C_Aq"]==0 | mut.cn.dat[, "SC_Aq_a"] == 0
   derived.LOH.ix = mut.cn.dat[, "C_Aq"]==0 | mut.cn.dat[, "SC_Aq_d"] == 0

   if( any(ancestral.LOH.ix) & mut.w[["GL"]] > 0) 
   {
#  TODO: still need to adapt this
#      GL_PR = rowSums(mut.ev.mat[LOH.ix, c("alt_minor", "alt_major", "hom_alt")])
#      post_Prs[LOH.ix, "Pr_GL_som_HZ_alt"] <- mut.ev.mat[LOH.ix, "alt_major"] / GL_PR
#      post_Prs[LOH.ix, "Pr_GL_som_HZ_ref"] <- mut.ev.mat[LOH.ix, "alt_minor"] / GL_PR
   }

   post_Prs[,"Pr_somatic_clonal"] = mut.ev.mat[,1]
   post_Prs[,"Pr_ge2"] = mut.ev.mat[,1] * rowSums(som_mut_Q_tab[,-1, drop=FALSE])

   post_Prs[,"Pr_subclonal"] = rowSums(mut.ev.mat[,c(2,3), drop=FALSE] )
   post_Prs[,"Pr_germline"] = mut.ev.mat[,4]

## som_mut_Q_tab refers to ancestral mut multiplicity q_a (= 0 for H3)
## H1: Ancestral homozygous iff  q_a == SC_Aq_a & C_Aq == 0  
##     Derived homozygous iff total_qs == C_Aq + SC_Aq_d == q_a == q_d
## H2: Ancestral homozygous iff  q_a == SC_Aq_a & C_Aq == 0 
##     Derived homozygous iff ancestral homozygous
## H3: SSNV is derived
##     Derived homozygous iff total_qs == q_d 
## H4: SSNV and SCNA in branched subpops;  
##     Derived homozygous iff total_qc == 1

   Pr_clonal = post_Prs[,"Pr_somatic_clonal"]
   Pr_subclonal = post_Prs[,"Pr_subclonal"]
## calculation for ancestral homozygosity
   Pr_q_a = som_mut_Q_tab
   r = (1:nrow(Pr_q_a)-1)
   total_qc = mut.cn.dat[, "total_qc"]
   total_qs = mut.cn.dat[, "total_qs"]

   kQ = dim(H.123_Pr_q)[3]
   N_SSNV = dim(H.123_Pr_q)[2]

#   H1_Pr_q_a = H.123_Pr_q[1,,]
#   H2_Pr_q_a = H.123_Pr_q[2,,]
#   H3_Pr_q_d = H.123_Pr_q[3,,]

   Pr_H1_anc_hom = rep(NA, N_SSNV)
   Pr_H2_anc_hom = rep(NA, N_SSNV)
   Pr_H1_der_hom = rep(NA, N_SSNV)
   Pr_H3_der_hom = rep(NA, N_SSNV)
   Pr_H4_der_hom = rep(0, N_SSNV)

## calculation for ancestral homozygosity
   for( i in 1:N_SSNV ) 
   {
      if( total_qc[i] == 0 ) { next }
      if( total_qc[i] < kQ ) { 
         Pr_H1_anc_hom[i] = H.123_Pr_q[1, i, total_qc[i] ] * H.ev[i,1]
         Pr_H2_anc_hom[i] = H.123_Pr_q[2, i, total_qc[i] ] * H.ev[i,2]
      } 
      else {
         Pr_H1_anc_hom[i] = 0
         Pr_H2_anc_hom[i] = 0
      }
   }
   
## calculation for derived homozygosity
   Pr_H2_der_hom = Pr_H2_anc_hom

   for( i in 1:N_SSNV ) 
   {
      if( total_qs[i] == 0 ) { next }

      if( total_qc[i] < kQ ) { 
         Pr_H1_der_hom[i] = H.123_Pr_q[1, i, total_qs[i] ] * H.ev[i,1]
         Pr_H3_der_hom[i] = H.123_Pr_q[3, i, total_qs[i] ] * H.ev[i,3]
      }
      else {
         Pr_H1_der_hom[i] = 0
         Pr_H3_der_hom[i] = 0
      }

      if( total_qc[i] == 1 ) { Pr_H4_der_hom[i] = H.ev[i,4] }
   }
   
   post_Prs[,"Pr_wt0"] = (Pr_H1_anc_hom + Pr_H2_anc_hom) * Pr_clonal

   post_Prs[,"Pr_subclonal_wt0"] = Pr_H1_der_hom + Pr_H3_der_hom + Pr_H4_der_hom + (Pr_H1_anc_hom + Pr_H2_der_hom * (1-Pr_clonal))

#   ix.01 = total_qc == 1 & total_qs == 0
#   post_Prs[ix.01, "Pr_subclonal_wt0"] = post_Prs[ix.01, "Pr_subclonal"]

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


calc_SSNV_on_subclonal_SCNA_CCF_dens = function( alpha, mut.cn.dat, allelic_subclonal_scna_tab, scna_log_ccf_dens, SSNV_model )
{
  get_phi_mat = function(grid, qc, qs, alpha) 
  {
    mat = matrix(NA, nrow=length(grid), ncol=length(grid) )
    f3 = grid
    for( i in 1:length(grid) )
    {
       f2 = grid[i]
       mat[i,] = (alpha / (2*(1-alpha) + alpha*qc*(1-f3) + alpha*qs*f3))
    }
    return(mat)
  }

  get_H3_phi_mat = function(grid, qc, qs, alpha) 
  {
    mat = matrix(NA, nrow=length(grid), ncol=length(grid) )
    f3 = grid
    for( i in 1:length(grid) )
    {
       f2 = grid[i]
       mat[i,] = (alpha / (2*(1-alpha) + alpha*qc*(1-f2-f3) + alpha*qs*(f2+f3)))
    }
    return(mat)
  }


  get_f2_plus_f3_dens_from_joint_f2_f3_dens = function( f2_f3_dens )
  {
    N_mut = dim(f2_f3_dens)[1]
    N_grid = dim(f2_f3_dens)[2]
    dd = matrix(NA, nrow=N_mut, ncol=N_grid )

    for( i in 1:N_mut )
    {
       dd[i,] = sum_marginalize_2D_sq_matrix( f2_f3_dens[i,,] )   [c(1:N_grid)]
    }

    return( dd )
  }

  get_f3_dens_from_joint_f2_f3_dens = function( f2_f3_dens )
  {
    dd = apply( f2_f3_dens, c(1,3), sum )
    return( dd )
  }

  get_f2_dens_from_joint_f2_f3_dens = function( f2_f3_dens )
  {
    dd = apply( f2_f3_dens, c(1,2), sum )
    return( dd )
  }


  N_grid = length(SSNV_model[["ccf_grid"]])
  H.123.ssnv.ccf.dens = array(NA, dim=c( 4, nrow(mut.cn.dat), N_grid) )
  H.123.ssnv.log.Z = array(NA, dim=c(nrow(mut.cn.dat), 4) ) 

  f2_grid = matrix(SSNV_model[["ccf_grid"]], nrow=N_grid, ncol=N_grid, byrow=FALSE) 
  f3_grid = matrix(SSNV_model[["ccf_grid"]], nrow=N_grid, ncol=N_grid, byrow=TRUE )

## put -Inf on CCF=1 (clonal)
  log_f2_prior = matrix( rep( log(1/(N_grid-1)), N_grid ), nrow=N_grid, ncol=N_grid, byrow=FALSE )
  log_f2_prior[N_grid,] = -Inf
 
  res = t(upper.tri(f2_grid, diag=TRUE))
  mask.ix = !res[ rev(c(1:nrow(res))), ]
  mask_mat = matrix(0, nrow=N_grid, ncol=N_grid )  ## renormalize so f2+f3 <= 1
  mask_mat[mask.ix] = -Inf

  for( i in 1:nrow(allelic_subclonal_scna_tab) )
  {
    SC.ix = allelic_subclonal_scna_tab[i,"SC.ix"]
    A1.ix = allelic_subclonal_scna_tab[i,"A1.ix"]
    A2.ix = allelic_subclonal_scna_tab[i,"A2.ix"]

    seg_SCNA_log_ccf_dens = scna_log_ccf_dens[SC.ix, , drop=TRUE] # prior on f3 for H1 and H2
    mut_ix = mut.cn.dat[,"A1.ix"] == A1.ix & mut.cn.dat[,"A2.ix"] == A2.ix 

## Joint prior for H1, H2, H4 - Pr(f3) is SCNA CCF.   Pr(f2) is unif
    log_f3_prior = matrix( seg_SCNA_log_ccf_dens, nrow=N_grid, ncol=N_grid, byrow=TRUE)
    joint_f2_f3_log_prior_mat = log_f2_prior + log_f3_prior + mask_mat
    joint_f2_f3_log_prior_mat = joint_f2_f3_log_prior_mat - LogAdd( as.vector(joint_f2_f3_log_prior_mat))

    alt = mut.cn.dat[mut_ix, "alt"]
    ref = mut.cn.dat[mut_ix, "ref"]
    tqc = as.integer(allelic_subclonal_scna_tab[i, "total_qc"])
    tqs = as.integer(allelic_subclonal_scna_tab[i, "total_qs"])

    phi_mat = get_phi_mat( SSNV_model[["ccf_grid"]], tqc, tqs, alpha )

## H1:
    if( tqs > 0 )
    {
       H1_AF_mat = (f2_grid + f3_grid) * phi_mat
       res = joint_SSNV_SCNA_grid_density(alt, ref, H1_AF_mat, joint_f2_f3_log_prior_mat, SSNV_model )
       H.123.ssnv.ccf.dens[1, mut_ix, ] = get_f2_plus_f3_dens_from_joint_f2_f3_dens( res[["f2_f3_dens"]] )
       H.123.ssnv.log.Z[mut_ix, 1] = res[["SSNV_log_Z"]]
    }
    else
    {
       H.123.ssnv.ccf.dens[1, mut_ix,] = 0
       H.123.ssnv.log.Z[mut_ix, 1] = -Inf
    }

## H2:
    q_a = 1  # by assumption
    if( tqc > tqs ) { q_d = 0 } # loss
    if( tqc < tqs ) { q_d = 2 } # gain
    tmp_mat = q_a * f2_grid + q_d * f3_grid
    z.ix = which(tmp_mat==0)
    H2_AF_mat = tmp_mat * phi_mat
    H2_AF_mat[z.ix] = 0  ## get around 0 * Inf = NaN  errors
    res = joint_SSNV_SCNA_grid_density(alt, ref, H2_AF_mat, joint_f2_f3_log_prior_mat, SSNV_model )
    H.123.ssnv.ccf.dens[2, mut_ix,] = get_f2_plus_f3_dens_from_joint_f2_f3_dens( res[["f2_f3_dens"]] )
    H.123.ssnv.log.Z[mut_ix, 2] = res[["SSNV_log_Z"]]

## H3:
 ## prior on f2+f3 = SCNA CCF_dens
    if( tqs > 0 )
    {
       H3_phi_mat = get_H3_phi_mat( SSNV_model[["ccf_grid"]], tqc, tqs, alpha )
       H3_joint_f2_f3_log_prior_mat = sum_prior_on_2D_sq_matrix( seg_SCNA_log_ccf_dens )
       H3_AF_mat = f3_grid * H3_phi_mat
       res = joint_SSNV_SCNA_grid_density(alt, ref, H3_AF_mat, H3_joint_f2_f3_log_prior_mat, SSNV_model )
       H.123.ssnv.ccf.dens[3, mut_ix,] = get_f3_dens_from_joint_f2_f3_dens( res[["f2_f3_dens"]] )
       H.123.ssnv.log.Z[mut_ix, 3] = res[["SSNV_log_Z"]]
    } 
    else
    {
       H.123.ssnv.ccf.dens[3, mut_ix,] = 0
       H.123.ssnv.log.Z[mut_ix, 3] = -Inf
    }

    
## H4: prior on f3 = SCNA CCF dens
#       q_a = 0; q_d=1
    H4_AF_mat = f2_grid * phi_mat
    res = joint_SSNV_SCNA_grid_density(alt, ref, H4_AF_mat, joint_f2_f3_log_prior_mat, SSNV_model )
    H.123.ssnv.ccf.dens[4, mut_ix, ] = get_f2_dens_from_joint_f2_f3_dens( res[["f2_f3_dens"]] )
    H.123.ssnv.log.Z[mut_ix, 4] = res[["SSNV_log_Z"]]

  }

  return( list( "H.123.ssnv.ccf.dens"= H.123.ssnv.ccf.dens, "H.123.ssnv.log.Z"=H.123.ssnv.log.Z) )
}
