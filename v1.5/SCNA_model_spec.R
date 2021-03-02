## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


## Totally generic
SCNA_model_eval_prior = function( SCNA_model )
{
  norm_LL = rep( NA, length(SCNA_model[["Theta"]]) )
## add Gaussian priors over Theta (transformed)
  for( i in 1:length(SCNA_model[["Theta"]]) )
  {
     norm_LL[i] =                dnorm( SCNA_model[["Theta"]][i],  
		 		 SCNA_model[["pi.Theta"]][[i]]["mu"],	
		  		 SCNA_model[["pi.Theta"]][[i]]["sigma"], log=TRUE )
  }

## and Dirichlet priors on mixture weight marginal likelihood
  dir_LL = dDirichlet( SCNA_model[["mix.w"]], SCNA_model[["pi.mix.w"]], log=TRUE )

  LL = sum(norm_LL) + dir_LL

  return(LL)
}


SCNA_model_init = function(SCNA_model)
{
  for( i in 1:length(SCNA_model[["Theta"]]) )
  {
    SCNA_model[["Theta"]][i] = SCNA_model[["pi.Theta"]][[i]]["mu"]
  }

  SCNA_model[["Theta"]]["sigma.A"] = SCNA_model[["provisional.sigma.A"]]

  pc = SCNA_model[["pi.mix.w"]]
  ## Dir mode
  SCNA_model[["mix.w"]] = (pc - 1) / (sum(pc) - length(pc)) 

  SCNA_model[["Theta"]]["theta.Q"] = log(1000) ## approx uniform on CN comb

  return(SCNA_model)
}


grid_compose_LL = function( LL, log_pow_tab, log.ccf.post, dx )
{
   nc = ncol(log.ccf.post)
   nr = nrow(log.ccf.post)

   r2 = matrix(LL, ncol=nc, nrow=nr, byrow=TRUE)
#     r3 = r2 + log_pow_tab + log.ccf.post + log(dx)
   r3 = r2 + log.ccf.post  #  + log(dx)
   res_LL <- LogAdd(r3) 
   res_LL[ is.nan(res_LL) ] = -Inf  # all -Inf in LogAdd(r3) ...

   return(res_LL)
}

SCNA_model_setup = function( SCNA.argv, verbose=FALSE )
{
  SCNA_model = list()

 # Hyper-parameters
  SCNA_model[["pi.mix.w"]] = c( 50, 5, 2, 25 )   ## Order matters!
  names(SCNA_model[["pi.mix.w"]]) =c("clonal", "amp", "exp", "unif")

  SCNA_model[["pi.Theta"]] = list( "lambda"=c("mu"=log(100), "sigma"=2),
				   "theta.Q"=c("mu"=log(1), "sigma"=1),
				   "sigma.A"=c("mu"=NA, "sigma"=NA) )
 
## platform specific
  if( SCNA.argv[["copy_num_type"]] == "total" )
  {
     SCNA_model[["theta.0"]] = 1e-3
     SCNA_model[["pi.Theta"]][["sigma.A"]]["mu"] = log(0.05) 
     SCNA_model[["pi.Theta"]][["sigma.A"]]["sigma"] = 0.25

     SCNA_model[["Theta_lower"]] = log(c(10, 0.5, 0.001 ))
     SCNA_model[["Theta_upper"]] = log(c(125, 10, 0.1))

     SCNA_model[["provisional.sigma.A"]] = log(0.2)  ## std dev. of CAPSEG probes
  }
  if( SCNA.argv[["copy_num_type"]] == "allelic" )
  {
     SCNA_model[["theta.0"]] = NA
     SCNA_model[["pi.Theta"]][["sigma.A"]]["mu"] = log(1)
     SCNA_model[["pi.Theta"]][["sigma.A"]]["sigma"] = 0.3

     SCNA_model[["Theta_lower"]] = log(c(25, 0.5, 0.8 ))
     SCNA_model[["Theta_upper"]] = log(c(125, 10, 1))
     SCNA_model[["provisional.sigma.A"]] = log(1)
  }

  SCNA_model[["sigma.h"]] = SCNA.argv[["sigma.h"]]

  SCNA_model[["Theta"]] = rep(NA, 3)
  names(SCNA_model[["Theta"]]) = names(SCNA_model[["pi.Theta"]])
  SCNA_model[["Theta_sym"]] = c("%", "*", "~")
  SCNA_model[["Theta_tol"]] = c(1e-2, 1e-2, 1e-3)

 ## Other stuff
  SCNA_model[["ccf_grid"]] = seq( 0, 1, by=0.01)
  NGRID = length(SCNA_model[["ccf_grid"]])
  SCNA_model[["clonal_CCF_bins"]] = c( c(1:5), c((NGRID-5):NGRID))
  SCNA_model[["kQ"]] = 8 
  SCNA_model[["kTauDom"]] = c(SCNA.argv[["min.ploidy"]], SCNA.argv[["max.ploidy"]])
  SCNA_model[["unif_CCF_dom"]] = c(0, 1)
#  SCNA_model[["Pr_sub_zero_CN"]] = 2e-4
  SCNA_model[["amp_seg_Pr_M_b"]] = c(-1, 0)
  SCNA_model[["neg_seg_Pr_M_b"]] = c(-1, 0)
  SCNA_model[["Pr_bi.allelic.outlier"]] = 1

  SCNA_model[["provisional_clonal_Pr_threshold"]] = 0.5   ## > than this discluded from DP clustering
  SCNA_model[["seg_sem_thresh"]] = 0.1   ## > than this discluded from DP clustering

  SCNA_model[["opttol"]] = 0.001

# for provisional sweep
  SCNA_model[["kDom1"]] <- c(0, 1)
  SCNA_model[["kDom2"]] <- log(c(0.04, 1.05))


  SCNA_model[["N_DP_iter"]] = 100
  SCNA_model[["N_DP_burn"]] = 50
  SCNA_model[["Pi_k"]] = list("r"=10, "mu" = 3)

 # Initialization for povisional search on alpha, tau
  SCNA_model = SCNA_model_init(SCNA_model)

  if( verbose )
  {
     print("Initialized SCNA model:")
     print_model(SCNA_model)
  }

  return(SCNA_model)
}


get_SCNA_model_fit_plot_strings = function( mode.info )
{
  msg1 <- substitute( paste( purity(hat(alpha)) == x1, ", ",
			   "ploidy: ", hat(tau) == x2, ", ",
                           hat(tau)[g] == x3,  sep = ""),
                     list(
                          x1 = round(mode.info["alpha"], 2),
                          x2 = round(mode.info["tau"], 2),
                          x3 = round(mode.info["genome mass"], 2)))
  

#print(msg1)


  msg2 <- substitute(paste( sigma[H] == x1, ", ",
                            hat(sigma[A]) == x2, ", ",
   			    hat(theta[Z]) == x3, ", ", 
			    hat(theta[Q]) == x4, ", ",
			    hat(lambda) == x5, sep = ""), 
                     list( x1 = round(mode.info["sigma.h.hat"], 3),
                           x2 = round(mode.info["sigma.A.hat"], 3),
			   x3 = round(mode.info["theta.z.hat"], 2),
                           x4 = round(mode.info["theta.Q.hat"], 3),
                           x5 = round(mode.info["lambda.hat"], 2)))
#                           x4 = round(100 * mode.info["frac.het"], 2)))
  
#print(msg2)

  msg3 <- substitute(paste("SCNAs" == x6, ", ", "Kar" == x7, ", ", "SSNVs" == x8, 
                           sep = ""),
                     list(x6 = round(mode.info["SCNA_LL"], 2),
                          x7 = round(mode.info["Kar_LL"], 2),
                          x8 = round(mode.info["SSNV_LL"], 2)))

  return( list(msg1, msg2, msg3))
}

