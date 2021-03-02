## Full analysis of tCN CCF on subclonal SCNA is not implemented for tCN yet - model all SSNVs as on closest clonal SCNA


total_get_muts_nearest_clonal_scna <- function(mut.cn.dat, seg.q.tab, Q) {
  if (!all(c("mut_seg_ix") %in% colnames(mut.cn.dat)) ) {
    stop("wrong colnames")
  }
  
  muts.p.q <- array(NA, dim=c( nrow(mut.cn.dat), Q))
  
  for (i in seq_len(nrow(mut.cn.dat))) {
    muts.p.q[i, ] <- seg.q.tab[mut.cn.dat[i,"mut_seg_ix"], ]   
  }
  
  muts.q.hat <- apply(muts.p.q, 1, which.max ) - 1

  mut.cn.tab <- matrix( muts.q.hat, ncol=1 )
  colnames(mut.cn.tab) = "q_hat"
  
  return(mut.cn.tab) 
}



total_get_subclonal_scna_mut_ix = function(mut.cn.dat, subclonal_scna_tab) 
{
  wix1 = which(subclonal_scna_tab[, "subclonal_ix"] == 1)
  wix2 = which( !is.na(subclonal_scna_tab[, "qs"]) )

#  ix = (mut.cn.dat[, "A1.ix"] %in% wix1 | mut.cn.dat[, "A2.ix"] %in% wix1 ) &
#       (mut.cn.dat[, "A1.ix"] %in% wix2 & mut.cn.dat[, "A2.ix"] %in% wix2 ) 

  return(rep(0,nrow(mut.cn.dat)))
  
#  return(ix)
}




get_total_subclonal_scna_tab = function(mut.cn.dat, subclonal_scna_tab, seg.q.tab ) 
{
  # loop over HSCR segs
#  a_ix_keys = paste(mut.cn.dat[, "A1.ix"], mut.cn.dat[, "A2.ix"], sep="__" )
#  u_keys= unique(a_ix_keys)

  u_keys= c()

  n_keys = length(u_keys)
  mut_n_keys = rep(NA, n_keys)

  tabcols = c("SC.ix", "total_qc", "total_qs")
  subclonal_scna_tab = matrix( NA, nrow=length(u_keys), ncol=length(tabcols) )
  colnames(subclonal_scna_tab) = tabcols

  if( n_keys == 0 )  ## No muts on subclonal SCNAs!
  {
    subclonal.mut.tab = as.data.frame(matrix(NA, nrow=0, ncol=ncol(subclonal_scna_tab) ) )
   colnames(subclonal.mut.tab) = colnames(subclonal_scna_tab)
#   subclonal.mut.tab = subclonal.mut.tab[ , setdiff(colnames(subclonal.mut.tab), c("A1.ix","A2.ix"))]
   return( list("subclonal_scna_tab" = subclonal_scna_tab, "subclonal.mut.tab"=subclonal.mut.tab) )
  }
  
  stop( "This should not happen!")
#  return( list("subclonal_scna_tab" = subclonal_scna_tab, "subclonal.mut.tab"=subclonal.mut.tab) )
}
