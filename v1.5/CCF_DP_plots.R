
sample_trans_ccf_plot = function(mut_cols, grid, ccf_dens, sid, y_max) {
  x_max = length(grid)
  plot(0, type='n', bty='n', xaxt="n", main="", ylim=c(0, y_max), xlim=c(1, x_max), 
       xlab="Fraction of cancer cells with alteration", ylab="Density")
  axis(side=1, at=seq(0, x_max - 1, length=11), labels=seq(0, 1, by=0.1))
  title(main=sid, cex=par("cex"))
  
  trans = 1000 / nrow(ccf_dens) 
  if (trans < 20) { 
    trans = 20 
  }
  if (trans > 255) { 
    trans = 255 
  }
  
  for(i in seq_len(nrow(ccf_dens))) {
    cr = col2rgb(mut_cols[i])
    use_col = rgb(cr[1, 1], cr[2, 1], cr[3, 1], trans, maxColorValue=255)
    
    xx = c(seq_len(ncol(ccf_dens)), rev(seq_len(ncol(ccf_dens))))
    yy = c(rep(0, ncol(ccf_dens)), rev(ccf_dens[i,]))
    polygon( xx, yy, border=FALSE, col=use_col, add=TRUE)
  }
}


DP_mut_posterior_plots = function( alt_tab, loc_res, GRID, CCF_dens, SID, remove_unknown=FALSE )
{   
   if( remove_unknown )
   {
      nix = alt_tab[,"locus"] %in% c("", "Unknown")
      alt_tab = alt_tab[ !nix, ]
      loc_res = loc_res[ !nix, ]
      CCF_dens = CCF_dens[ !nix, ]
   }

## sort mutations by 1st order of posterior CCF dist
   E_loc = colSums( GRID * t(loc_res) )
   loc_res = loc_res[ order(E_loc), ] 
   CCF_dens = CCF_dens[ order(E_loc), ]
   alt_tab = alt_tab[ order(E_loc), ]


   dmar = par("mar")
   nmar=dmar
   nmar[c(3)] = 2  # small top margin
   nmar[c(1)] = 0  # no bot margin
   par(mar=nmar)
   plot.new()
   title(main=SID, cex=par("cex"))
   
   nmar=dmar
   nmar[c(1,3)] = 0  # no bot or top margin
   par(mar=nmar)
   
   #plot( 0, type="n", bty="n", main="", xlab="Simulated CCF", ylab="Density", xlim=c(0,1), ylim=c(0, ymax), las=1 )
   
   for( i in 1:nrow(loc_res) )
   {
      mut_col = ifelse( alt_tab[i,"type"] == "SCNA", "cyan4", "maroon" )
#      font = ifelse( alt_tab[i,"type"] == "SCNA", 1, 3 )
      font = 1
   
      plot( 0, type="n", bty="n", axes=FALSE, main="", ylab="", xlab="", xlim=c(1,length(GRID)), ylim=c(0,0.25) )

      barplot( loc_res[i,], space=0, axisnames=FALSE, yaxt="n", border=FALSE, col=mut_col, add=TRUE )
   
      lines( c(1:length(GRID)), CCF_dens[i,], col="grey15", lwd=2 )
   
      if(  alt_tab[i,"type"] == "SCNA" )
      {
         mut_lab = alt_tab[i,"locus"] 
      }
      if(  alt_tab[i,"type"] == "SNV" )
      {
         var_info = ifelse( alt_tab[i,"Variant_Classification"] %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "In_Frame_Ins"), as.character(alt_tab[i,"Protein_Change"]),  as.character(alt_tab[i,"Variant_Classification"]) )
         cov_1d = paste( alt_tab[i,"alt"], "/", alt_tab[i, "coverage"], sep="")

         mut_lab_1 = paste( alt_tab[i,"locus"], ", ", var_info, sep="" )
         mut_lab_2 = paste( alt_tab[i,"mut_class"], ", (", cov_1d, ")", sep="" ) 

         if( is.na(alt_tab[i,"H1"]) ) {
            l1 = paste( alt_tab[i,"HS_q_hat_1"], ",", alt_tab[i,"HS_q_hat_2"], sep="" )
         }
         else {
            v = alt_tab[i,c("SC_Aq_a", "SC_Aq_d", "C_Aq")]
            l1 = paste( v[1], ">", v[2], ",", v[3], sep="")
         }
         mut_lab_3 = paste( "CN: ", l1, sep="")
      }

      mut_lab = paste( mut_lab_1, mut_lab_2, mut_lab_3, sep="; ")
      mtext( text=mut_lab, side=3, line=-1, adj=0, cex=par("cex"), font=font )
   }
   
   nmar[1] = 3.5
   nmar[3] = 0 
   par(mar=nmar)
   dmgp=par("mgp")
   nmgp=dmgp
   nmgp[1]=2  ## move axis title closer in
   par(mgp=nmgp)
   
   plot( 0, type='n', bty='n', axes=FALSE, main="", ylab="", xlim=c(0,1), xlab="Fraction of cancer cells with alteration" )
   axis( side=1, at=seq(0,1,by=0.1) )
}


DP_post_gamma_plot = function( DP_res, N_burn )
{
   DP_post = DP_res$DP_post
   DP_prior = DP_res$DP_prior
   Pi_gamma = DP_prior$Pi_gamma

   N_mut = dim(DP_post$assign)[1]
   N_iter = length(DP_post$c_fpost)
   est_gamma = rep( NA, N_iter-N_burn+1 )
   XMAX = 5 
   GRID = seq(1e-5, XMAX, length=500)
   gamma_post = array( NA, dim=c(N_iter-N_burn+1, length(GRID)) )

   Pi_gamma = DP_prior$Pi_gamma
   a = Pi_gamma$a  
   b = Pi_gamma$b

   for( i in N_burn:N_iter )
   {
      est_gamma[(i-N_burn+1)] = DP_post$gamma[i]
      k = DP_post$K[i]
      eta = DP_post$eta[i]
      m1 = dgamma( GRID, a + k, b - log(eta) ) 
      m2 = dgamma( GRID, a + k - 1, b - log(eta) )
      D = N_mut * (b - log(eta))
      w = (a+k-1) / (D+a+k)
      gamma_post[(i-N_burn+1),] = w * m1 + (1-w) * m2
   }

   gamma_sampled = hist( est_gamma, breaks=40, plot=FALSE)$density

   gamma_mat = rbind( dgamma(GRID, a,b), colMeans(gamma_post)  )

   YMAX = max( c(gamma_mat, gamma_sampled) )
   XMAX = GRID[  max( which( colSums(gamma_mat) > 1e-3)) ]

   plot( 0, type="n", ylim=c(0,YMAX), xlim=c(0.01,XMAX), lty=2, bty="n", xlab="DP concentration (gamma)", ylab="Density", main="" )
   legend( x='topright', legend=c("P(gamma)", "P(gamma | D)"), lty=c(1), col=c(1,2), bty="n" )

   hist( est_gamma, breaks=40, add=TRUE, freq=FALSE, col="grey", border=FALSE )

   lines( GRID, gamma_mat[1,], col=1 )
   lines( GRID, gamma_mat[2,], col=2 )
}






## express prior and posterior as distribution over number of DP components (K).
DP_post_K_plot= function( DP_res, N_burn, GRID )
{
   DP_post = DP_res$DP_post
   DP_prior = DP_res$DP_prior
   k_0_map = DP_prior$k_0_map
   Pi_gamma = DP_prior$Pi_gamma

   N_mut = dim(DP_post$assign)[1]
   a = Pi_gamma$a  
   b = Pi_gamma$b

   N_iter = length(DP_post$c_fpost)
   est_K = rep( NA, N_iter-N_burn+1 )

   for( i in N_burn:N_iter )
   {
      est_K[(i-N_burn+1)] = DP_post$K[i]
   }

   est_k_post = rep( 0, N_iter-N_burn+1 )
   for( i in N_burn:N_iter )
   {
      k = DP_post$K[i]
      eta = DP_post$eta[i]
      gamma = DP_post$gamma[i]
      m1 = dgamma( GRID, a + k, b - log(eta) ) 
      m2 = dgamma( GRID, a + k - 1, b - log(eta) )
      D = N_mut * (b - log(eta))
      w = (a+k-1) / (D+a+k)
 
      MCMC_ix = i-N_burn+1

      GRID_ix = which.min( (k_0_map[,2] - gamma)^2 )
      est_k_post[MCMC_ix] = k_0_map[GRID_ix, 1]
   }
   k_post = hist( est_k_post, breaks=c(1:(N_mut+1))-0.5, plot=FALSE)$density

   max1 = max( which( k_post > 1e-3 ) )
   max2 = max( k_0_map[,1][ DP_prior$k_prior > 1e-3 ] )
   max3 = max( est_k_post )
   XMAX= max(max1, max2, max3)

# approximated prior
   DG = dgamma( k_0_map[,2],  DP_prior$Pi_gamma$a, DP_prior$Pi_gamma$b ) 
   k_prob = DG
   k_prob = k_prob/sum(k_prob) * c( 1, 1/diff(c(k_0_map[,1])) )

   YLIM = max( c(DP_prior$k_prior, k_prob, k_post) )

   plot( 0, type='n', bty="n",  xlab="Number of clusters (k)", ylab="Density", main="", ylim=c(0,YLIM), xlim=c(0, XMAX) )

 # sampled k counts
   hist( est_K, breaks=c(1:(N_mut+1))-0.5, add=TRUE, freq=FALSE, col="grey", border=FALSE )

 # N-binom prior
   lines( c(1:N_mut), DP_prior$k_prior, col=1 )

# approximated prior
   lines( k_0_map[,1], k_prob, col=1, lty=2 )

 # predicted posterior
   lines( c(1:length(k_post)), k_post, col=2 )


   legend( x='topright', legend=c("P(k)", "P(k) approx", "P(k | D)", "sampled k values"), col=c(1,1,2,"grey"), bty="n", lty=c(1,2,1,1), lwd=2 )

   mtext( text=paste("N=", N_mut,sep=""), cex=par("cex"), adj=0)

}



