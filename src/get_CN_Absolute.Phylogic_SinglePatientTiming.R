library(plyr)
library(optparse)
option.list <- list(
	make_option("--rdata_fn", dest="Called_Absoute"),
	make_option("--pair_id", dest="pair_id"),
	make_option("--abs_lib_dir", dest="abs_lib_dir"),
	make_option("--pcawg", dest="pcawg", default="FALSE"))

opt <- parse_args(OptionParser(option_list=option.list))
print(opt)

Called_Absolute= opt[["Called_Absoute"]]
CGA_DIR_ABS= opt[["abs_lib_dir"]]
pair_id= opt[["pair_id"]]
PCAWG= opt[["pcawg"]]

Out_Name=file.path(paste(pair_id,"tsv",sep="."))

print( paste("sourcing files in ", CGA_DIR_ABS, sep=""))

rr = dir(CGA_DIR_ABS,full.names=TRUE, pattern="*.R$")
for( i in 1:length(rr) ) {
  source(rr[i])}
load(Called_Absolute)
input_purity=seg.obj$mode.res$mode.tab[[1]]

allele.segs= get_hom_pairs_segtab(seg.obj)

ExtractSampleObs <<- AllelicExtractSampleObs
GetCopyRatioComb <<- AllelicGetCopyRatioComb
SCNA_model = seg.obj[["mode.res"]][["mode_SCNA_models"]][[1]] 
tree_clust = resort_tree_clusters( SCNA_model, SCNA_model[["seg_CCF_DP"]][["tree_clust"]] )

SCNA_model[["seg_CCF_DP"]][["tree_clust"]] = tree_clust
SCNA_model[["seg_CCF_DP"]][["seg_clust_tab"]] = get_seg_clust_tab( SCNA_model )
mut_cols = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["assign"]]

mode.tab <- seg.obj[["mode.res"]][["mode.tab"]]

res = get_b_and_delta( mode.tab[1, "alpha"], mode.tab[1, "tau"] )
delta = res$delta
b = res$b
segs_d0 = deconstruct_SCNAs( SCNA_model, seg.obj, allele.segs, b, delta )

AS.seg.ix = allele.segs[, c("seg.ix.1", "seg.ix.2")]
d0.allele.segs = allele.segs
d0.allele.segs[,"A1.Seg.CN"] =segs_d0[ AS.seg.ix[,1] ]
d0.allele.segs[,"A2.Seg.CN"] = segs_d0[ AS.seg.ix[,2] ]
 
d0.allele.segs["seg.ix.1"] = NULL
d0.allele.segs["seg.ix.2"] = NULL
d0.allele.segs["W"] = NULL
d0.allele.segs=rename(d0.allele.segs, c("Start.bp"="Start", 
"End.bp"="End","n_probes"="N_probes","length"="Length","A1.Seg.sem"="A1.Sigma","A2.Seg.sem"="A2.Sigma")) 

write.table(d0.allele.segs, file = Out_Name, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)
if(PCAWG=="TRUE"){
Out_Name=file.path(paste(pair_id,"segments.txt",sep="_"))
load(Called_Absolute)
input_purity=seg.obj$mode.res$mode.tab[[1]]

  ExtractSampleObs <<- AllelicExtractSampleObs
  GetCopyRatioComb <<- AllelicGetCopyRatioComb
  SCNA_model = seg.obj[["mode.res"]][["mode_SCNA_models"]][[1]]
  tree_clust = resort_tree_clusters( SCNA_model, SCNA_model[["seg_CCF_DP"]][["tree_clust"]] )

  SCNA_model[["seg_CCF_DP"]][["tree_clust"]] = tree_clust
  SCNA_model[["seg_CCF_DP"]][["seg_clust_tab"]] = get_seg_clust_tab( SCNA_model )
  mut_cols = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["assign"]]

  mode.tab <- seg.obj[["mode.res"]][["mode.tab"]]

  res = get_b_and_delta( mode.tab[1, "alpha"], mode.tab[1, "tau"] )
  delta = res$delta
  b = res$b
  segs_d0 = deconstruct_SCNAs( SCNA_model, seg.obj, allele.segs, b, delta )


  sc_tab = seg.obj[["mode.res"]][["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]][1,,]

  nc = ncol(SCNA_model[["collapsed_CCF_dens"]])
  GRID = cumsum( c(0, rep(1/(nc-1),(nc-1))))
  clust_dens = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["CCF_dens"]]
  clust_dens_str = apply(clust_dens,1,paste,collapse=",")
  clonal_clust_dens_str = paste(as.character(c(0:89*0,1)), collapse=",")
  #clust_dens_str = paste(as.character(one), collapse=", ")
  ccf_hats=max.col(clust_dens)/100
  mut_cols = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["assign"]]

  sc_ix=mut_cols != which(clust_dens[,1]==1)
  sc_ix[is.na(sc_ix) & sc_tab[,"subclonal_ix"]==0]=0
  sc_ix[sc_ix==FALSE]=0
  sc_ix[sc_ix==TRUE]=1

  #sc_ix_noNA=sc_ix
  #sc_ix_noNA[is.na(sc_ix)]=-1
  #segs_d0[sc_ix_noNA==0]=round(segs_d0[sc_ix_noNA==0])

  AS.seg.ix = allele.segs[, c("seg.ix.1", "seg.ix.2")]
  d0.allele.segs = allele.segs
  d0.allele.segs[,"A1.Seg.CN"] = round(segs_d0[ AS.seg.ix[,1] ])
  d0.allele.segs[,"A2.Seg.CN"] = round(segs_d0[ AS.seg.ix[,2] ])

  d0.allele.segs["W"] = NULL
  d0.allele.segs=rename(d0.allele.segs, c("Start.bp"="Start", 
"End.bp"="End","n_probes"="N_probes","length"="Length","A1.Seg.sem"="A1.Sigma","A2.Seg.sem"="A2.Sigma"))

  d0.allele.segs["seg.ix.1"] = NULL
  d0.allele.segs["seg.ix.2"] = NULL
  d0.allele.segs["A1.Sigma"] = NULL
  d0.allele.segs["A2.Sigma"] = NULL

  A1.Subclonal.ix = sc_ix[AS.seg.ix[,1]]
  A2.Subclonal.ix = sc_ix[AS.seg.ix[,2]]
  A1.Subclonal.CN = sc_tab[,"qs"][AS.seg.ix[,1]]
  A2.Subclonal.CN = sc_tab[,"qs"][AS.seg.ix[,2]]
  A1.Clonal.CN = sc_tab[,"qc"][AS.seg.ix[,1]]
  A2.Clonal.CN = sc_tab[,"qc"][AS.seg.ix[,2]]

  A1.CCF_hat = ccf_hats[mut_cols[AS.seg.ix[,1]]]
  A2.CCF_hat = ccf_hats[mut_cols[AS.seg.ix[,2]]]
  A1.CCF_dens = clust_dens_str[mut_cols[AS.seg.ix[,1]]]
  A2.CCF_dens = clust_dens_str[mut_cols[AS.seg.ix[,2]]]


  A1.CCF_hat[is.na(A1.CCF_hat)]=sc_tab[,"CCF_hat"][AS.seg.ix[,1]][is.na(A1.CCF_hat)]
  A2.CCF_hat[is.na(A2.CCF_hat)]=sc_tab[,"CCF_hat"][AS.seg.ix[,2]][is.na(A2.CCF_hat)]

  d0.allele.segs[,"CCF"]=1
  d0.allele.segs[,"CCF_dens"]=clonal_clust_dens_str

  clonal.ix=intersect(which(A1.Subclonal.ix == 0),which(A2.Subclonal.ix == 0))
  a1.sub.ix=intersect(which(A1.Subclonal.ix != 0),which(A2.Subclonal.ix == 0))
  a2.sub.ix=intersect(which(A1.Subclonal.ix == 0),which(A2.Subclonal.ix != 0))
  all.sub.ix=intersect(which(A1.Subclonal.ix != 0),which(A2.Subclonal.ix != 0))
  All_Clonal_Seg=d0.allele.segs[intersect(which(A1.Subclonal.ix == 0),which(A2.Subclonal.ix == 0)),]
  A1_Sub_Clonal_Seg=d0.allele.segs[intersect(which(A1.Subclonal.ix != 0),which(A2.Subclonal.ix == 0)),]
  A2_Sub_Clonal_Seg=d0.allele.segs[intersect(which(A1.Subclonal.ix == 0),which(A2.Subclonal.ix != 0)),]
  Both_Sub_Clonal_Seg=d0.allele.segs[intersect(which(A1.Subclonal.ix != 0),which(A2.Subclonal.ix != 0)),]


  d0.allele.segs[,"Historically.Clonal"]=0

  if (length(a1.sub.ix) > 0){
    d0.allele.segs[a1.sub.ix,"A1.Seg.CN"]=A1.Clonal.CN[a1.sub.ix]
    d0.allele.segs[a1.sub.ix,"CCF"]=1#-A1.CCF_hat[a1.sub.ix]
    d0.allele.segs[a1.sub.ix,"CCF_dens"]=clonal_clust_dens_str
    d0.allele.segs[a1.sub.ix,"Historically.Clonal"]=1


    a1_segs=d0.allele.segs[a1.sub.ix,]
    a1_segs["A1.Seg.CN"]=A1.Subclonal.CN[a1.sub.ix]
    a1_segs["CCF"]=A1.CCF_hat[a1.sub.ix]
    a1_segs["CCF_dens"]=A1.CCF_dens[a1.sub.ix]
    a1_segs["Historically.Clonal"]=0

  }else{
  a1_segs=NULL}
  if (length(a2.sub.ix) > 0){
    d0.allele.segs[a2.sub.ix,"A2.Seg.CN"]=A2.Clonal.CN[a2.sub.ix]
    d0.allele.segs[a2.sub.ix,"CCF"]=1#-A2.CCF_hat[a2.sub.ix]
    d0.allele.segs[a2.sub.ix,"CCF_dens"]=clonal_clust_dens_str
    d0.allele.segs[a2.sub.ix,"Historically.Clonal"]=1

    a2_segs=d0.allele.segs[a2.sub.ix,]
    a2_segs["A2.Seg.CN"]=A2.Subclonal.CN[a2.sub.ix]
    a2_segs["CCF"]=A2.CCF_hat[a2.sub.ix]
    a2_segs["CCF_dens"]=A2.CCF_dens[a2.sub.ix]
    a2_segs["Historically.Clonal"]=0
  }else{

  a2_segs=NULL}
  if (length(all.sub.ix) > 0){

    d0.allele.segs[all.sub.ix,"A1.Seg.CN"]=A1.Clonal.CN[all.sub.ix]
    d0.allele.segs[all.sub.ix,"A2.Seg.CN"]=A2.Clonal.CN[all.sub.ix]
    matched= A1.CCF_hat[all.sub.ix] == A2.CCF_hat[all.sub.ix]

    d0.allele.segs[all.sub.ix[!matched],"CCF"]=1#-A1.CCF_hat[all.sub.ix[!matched]]-A2.CCF_hat[all.sub.ix[!matched]]
    d0.allele.segs[all.sub.ix[!matched],"CCF_dens"]=clonal_clust_dens_str
    d0.allele.segs[all.sub.ix[!matched],"Historically.Clonal"]=1

    d0.allele.segs[all.sub.ix[matched],"CCF"]=1#-A1.CCF_hat[all.sub.ix[matched]]
    d0.allele.segs[all.sub.ix[matched],"CCF_dens"]=clonal_clust_dens_str
    d0.allele.segs[all.sub.ix[matched],"Historically.Clonal"]=1
    if (length(all.sub.ix[!matched])> 0){
      all_sub_a1_segs=d0.allele.segs[all.sub.ix[!matched],]
      all_sub_a1_segs["A1.Seg.CN"]=A1.Subclonal.CN[all.sub.ix[!matched]]
      all_sub_a1_segs["CCF"]=A1.CCF_hat[all.sub.ix[!matched]]
      all_sub_a1_segs["CCF_dens"]=A1.CCF_dens[all.sub.ix[!matched]]
      all_sub_a1_segs["Historically.Clonal"]=0

      all_sub_a2_segs=d0.allele.segs[all.sub.ix[!matched],]
      all_sub_a2_segs["A2.Seg.CN"]=A2.Subclonal.CN[all.sub.ix[!matched]]
      all_sub_a2_segs["CCF"]=A2.CCF_hat[all.sub.ix[!matched]]
      all_sub_a2_segs["CCF_dens"]=A2.CCF_dens[all.sub.ix[!matched]]
      all_sub_a2_segs["Historically.Clonal"]=0
    }else{

      all_sub_a1_segs=NULL
      all_sub_a2_segs=NULL

    }
    if (length(all.sub.ix[matched])> 0){
      euqal_sub_seg=d0.allele.segs[all.sub.ix[matched],]
      euqal_sub_seg["A1.Seg.CN"]=A1.Subclonal.CN[all.sub.ix[matched]]
      euqal_sub_seg["A2.Seg.CN"]=A2.Subclonal.CN[all.sub.ix[matched]]
      euqal_sub_seg["CCF"]=A1.CCF_hat[all.sub.ix[matched]]
      euqal_sub_seg["CCF_dens"]=A1.CCF_dens[all.sub.ix[matched]]
      euqal_sub_seg["Historically.Clonal"]=0

    }else{

    euqal_sub_seg=NULL}}

  else{

    all_sub_a1_segs=NULL
    all_sub_a2_segs=NULL
    euqal_sub_seg=NULL
  }   

comb_segtab=rbind(d0.allele.segs,a1_segs,a2_segs,all_sub_a1_segs,all_sub_a2_segs,euqal_sub_seg)

comb_segtab=comb_segtab[with(comb_segtab, order(as.integer(comb_segtab[["Chromosome"]]), as.integer(comb_segtab[["Start"]]))), ]

comb_segtab[,"copy_number"]=comb_segtab[["A1.Seg.CN"]]+comb_segtab[["A2.Seg.CN"]]
comb_segtab["N_probes"]=NULL
comb_segtab["Length"]=NULL

comb_segtab=rename(comb_segtab, c("Chromosome"="chromosome","Start"="start", 
"End"="end","A1.Seg.CN"="minor_cn","A2.Seg.CN"="major_cn","CCF"="ccf","Historically.Clonal"="historically_clonal")) 
comb_segtab=comb_segtab[,c(1,2,3,9,4,5,6,8)]

write.table(comb_segtab, file = Out_Name, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)
}

