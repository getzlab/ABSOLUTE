## infer events from observed states

deconstruct_mode_results_SCNAs = function( seg.dat, mode.res, mut.cn.dat, verbose=verbose )
{
   allele.segs = get_hom_pairs_segtab( seg.dat )
   N = length( mode.res[["mode_SCNA_models"]] )

   for( i in 1:N )
   {
      b = mode.res[["mode.tab"]][i,"b"]
      delta = mode.res[["mode.tab"]][i,"delta"]
   
      mode.res[["mode_SCNA_models"]][[i]][["segs_d0"]] = deconstruct_SCNAs( mode.res[["mode_SCNA_models"]][[i]], seg.dat, allele.segs, b, delta )
   }

   return(mode.res)
}

deconstruct_SCNAs = function( SCNA_model, seg.dat, allele.segs, b, delta )
{ 
   CN_states = SCNA_model[["CN_states"]]
   seg.post_mix.w = SCNA_model[["seg.post_mix.w"]]

   obs = seg.dat[["obs.scna"]]
#   allele.segs = cbind( allele.segs, obs[["AS.seg.ix"]] )

   hom_seg.post_mix.w = array( NA, dim=c( 2, nrow(allele.segs), ncol(seg.post_mix.w) ) )
   hom_seg.post_mix.w[1,,] =  seg.post_mix.w[ allele.segs[, "seg.ix.1"], ]
   hom_seg.post_mix.w[2,,] =  seg.post_mix.w[ allele.segs[, "seg.ix.2"], ]

   tree_clust = SCNA_model[["seg_CCF_DP"]][["tree_clust"]]
   seg_clust_tab = SCNA_model[["seg_CCF_DP"]][["seg_clust_tab"]]
   n.ix = apply( SCNA_model[["seg.ix.tab"]][ , c("amp.ix", "neg.ix" ) ], 1, any )
# 1st cluster is clonal

# debug tree vs.  post
#   seg_modal_clusts = apply(seg_clust_tab[!n.ix,], 1, which.max)
#   same = sum(tree_clust[["assign"]][!n.ix] == seg_modal_clusts, na.rm=T)

## clonal segs collapse to CN_states[,"qc"]
## sc segs go to modal DP assignment X CN_states

   seg_modal_clust = rep( NA, nrow(seg_clust_tab) )
   if( any( !n.ix ) )
   {
      seg_modal_clust[!n.ix] = apply( seg_clust_tab[!n.ix,, drop=FALSE], 1, which.max )
   }

   clonal_clust_num = tree_clust[["CCF_order"]][1]

   clonal.ix = seg_modal_clust == clonal_clust_num 
   clonal.ix[is.na(clonal.ix)] = FALSE
   ng = ncol(SCNA_model[["seg_CCF_dens"]])
   clonal.bit = rep(FALSE, length(clonal.ix) )
   clonal.bit[clonal.ix] = SCNA_model[["seg_CCF_dens"]][clonal.ix,1] > SCNA_model[["seg_CCF_dens"]][clonal.ix,ng]

   subclonal.ix = !n.ix & !clonal.ix
   subclonal.ix[is.na(subclonal.ix)] = FALSE

   segs_d0 = rep(NA, length(clonal.ix) )
   segs_d0[clonal.ix & clonal.bit] = CN_states[ clonal.ix & clonal.bit , "qc"]
   segs_d0[clonal.ix & !clonal.bit] = CN_states[ clonal.ix & !clonal.bit , "qs"]

   s.qc = CN_states[subclonal.ix, "qc"]  
   s.qs = CN_states[subclonal.ix, "qs"]  

   s.CCF = tree_clust[["CCF_dens"]]  [ tree_clust[["assign"]][subclonal.ix], -1, drop=FALSE]   # remove clonal bin (1)
# compute expectation on subclonal CCFs using collapsed CCF grid 
   ccf_grid = SCNA_model[["ccf_grid"]] 
   sc_ccf_grid = ccf_grid[ - SCNA_model[["clonal_CCF_bins"]] ]
   s.E.CCF = rowSums( s.CCF * matrix( sc_ccf_grid, nrow=nrow(s.CCF), ncol=length(sc_ccf_grid), byrow=TRUE) )

   segs_d0[subclonal.ix] =  s.qc + s.E.CCF * (s.qs - s.qc)

   amp.ix = SCNA_model[["seg.ix.tab"]][,"amp.ix"]
   amp.CR = obs[["d.tx"]][amp.ix]
   neg.ix = SCNA_model[["seg.ix.tab"]][,"neg.ix"]
   segs_d0[neg.ix] = 0
   segs_d0[amp.ix] = (amp.CR-b) / delta ## rescale for PP

#   bix = which( seg_modal_clust != tree_clust[["assign"]] )

   return(segs_d0)


   hom_seg_clust_tab =  array( NA, dim=c( 2, nrow(allele.segs), ncol(seg_clust_tab) ) )
   hom_seg_clust_tab[1,,] = seg_clust_tab[ allele.segs[, "seg.ix.1"], ]
   hom_seg_clust_tab[2,,] = seg_clust_tab[ allele.segs[, "seg.ix.2"], ]


   hom_seg_event_tab =  matrix( NA, nrow=nrow(allele.segs), ncol=2 )
   hom_seg_event_tab[ clonal.ix, 1] = CN_states[ allele.segs[, "seg.ix.1"], "qc" ]
   hom_seg_event_tab[ clonal.ix, 2] = CN_states[ allele.segs[, "seg.ix.2"], "qc" ]

   chrs = unique(allele.segs[,"Chromosome"])
   for( i in seq_along(chrs))
   {
      chr.ix = allele.segs[,"Chromosome"] == chrs[i]

## indices into seg.mix tab
      A1.wix = allele.segs[chr.ix, "seg.ix.1"] 
      A2.wix = allele.segs[chr.ix, "seg.ix.2"]
   }
}


