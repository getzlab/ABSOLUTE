## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


SCNA_model_opt_cond_b_delta = function( obs, b, delta, SCNA_model, max.iter=25, verbose=FALSE )
{
   cur.loglik = Theta_LL_wrapper( SCNA_model[["Theta"]], obs, b, delta, SCNA_model )
   print( paste("Initial LL = ", round(cur.loglik,4), sep=""))
   print_model(SCNA_model)

   iter=1
   while( 1 )
   {
      SCNA_model[["Theta"]] = SCNA_model_Theta_sweep_coord_ascent( obs, b, delta, SCNA_model, verbose=verbose ) 

      SCNA_model[["mix.w"]] = SCNA_model_opt_mix.w( obs, b, delta, SCNA_model, verbose=FALSE) 

      loglik = SCNA_model_loglik( obs, b, delta, SCNA_model )
      cond <- abs(cur.loglik - loglik) / abs(cur.loglik)

      print(paste("loglik = ", round(loglik,4), sep=""))
      print(paste("cond = ", round(cond, 4), sep=""))
      print_model(SCNA_model)

      if(( cond < SCNA_model[["opttol"]]) || (iter >= max.iter)) {
         break
      }

      iter = iter+1
      cur.loglik = loglik
   }

   cat("\n")
   SCNA_model[["LL"]] = loglik 

   return(SCNA_model)
}




SCNA_Theta_1d_opt_step <- function( obs, b, delta, SCNA_model, parname, limits, symbol, opttol, verbose=FALSE) 
{
  LL <- function(par, parname) 
  {
    SCNA_model[["Theta"]][parname] <- par
    LL = SCNA_model_loglik( obs, b, delta, SCNA_model )

    if (!is.finite(LL)) 
    {
        print(paste(parname, ": Non-finite log-liklihood!", sep=""))
        stop()
    }
    else
    {
       if (verbose) { cat(symbol) }
    }

    return(-LL)
  }

  res <- optimize(LL, lower=limits[["lower"]], upper=limits[["upper"]],
                  tol=opttol, maximum=FALSE, parname=parname )

  return(res[["minimum"]])
}



SCNA_model_Theta_sweep_coord_ascent <- function( obs, b, delta, SCNA_model, verbose=FALSE) 
{
  symbols = SCNA_model[["Theta_sym"]]  
  tol = SCNA_model[["Theta_tol"]]
  lower = SCNA_model[["Theta_lower"]]
  upper = SCNA_model[["Theta_upper"]]
  Theta = SCNA_model[["Theta"]]

  for( i in 1:length(Theta) )
  {
    Theta[i] <- SCNA_Theta_1d_opt_step( obs, b, delta, SCNA_model, names(Theta)[i], list("lower"=lower[i], "upper"=upper[i]), symbols[i], tol[i], verbose=verbose)

    SCNA_model[["Theta"]][i] = Theta[i]
  }
  cat("\n")

  return( Theta )
}

SCNA_model_opt_mix.w <- function( obs, b, delta, SCNA_model, verbose=FALSE) 
{
  seg_mix_post = SCNA_model_seg_mix.w_post( obs, b, delta, SCNA_model )

  dc = colSums(seg_mix_post)
  pc = dc + SCNA_model[["pi.mix.w"]]

  ## Dir mode
  mode.w = (pc - 1) / (sum(pc) - length(pc)) 
  names(mode.w) = names( SCNA_model[["pi.mix.w"]])

  return(mode.w)
}




Theta_LL_wrapper = function( Par, obs, b, delta, SCNA_model )
{
   SCNA_model[["Theta"]] = Par
   LL = SCNA_model_loglik( obs, b, delta, SCNA_model )

   return(LL)
}


SCNA_model_Theta_Laplace_approx = function( obs, b, delta, SCNA_model, verbose=verbose )
{
## include b and delta in laplacian
  LL_wrapper = function( Par, obs, SCNA_model )
  {
     theta_names = names(SCNA_model[["Theta"]])

     SCNA_model[["Theta"]] = Par[theta_names]
     b = Par["b"]
     delta = exp( Par["delta"] )

     LL = SCNA_model_loglik( obs, b, delta, SCNA_model )
     return(LL)
  }

## Theta already optimized
   Theta_hat = SCNA_model[["Theta"]]

#   rr = GetAlphaAndTau(b, delta) 
#   alpha = rr[["alpha"]]
#   tau = rr[["tau"]]

   Par = c(Theta_hat, b, log(delta) )

   LL <- LL_wrapper( Par, obs, SCNA_model )
   hess.mat <- hessian( LL_wrapper, x=Par, method = "Richardson", obs = obs, SCNA_model=SCNA_model )

   logI = Laplace_approx( hess.mat, LL )

   return(logI)
}

post_Dir_log_norm_const = function( SCNA_model )
{
  seg_mix_post = SCNA_model_seg_post( obs, b, delta, SCNA_model )

  dc = colsums(seg_mix_post)
  alpha = dc + SCNA_model[["pi.mix.w"]]

  log_num = sum( lgamma(alpha) )
  log_den = lgamma(sum(alpha))

  return( log_num - log_den )
}


SCNA_model_Theta_tx = function( Theta )
{
   res = exp(unlist(Theta)) 
   names(res)=names(Theta)

   return(res)
}


print_model = function(SCNA_model)
{
  Theta=SCNA_model[["Theta"]]
  msg = paste( names(Theta), " = ", round(SCNA_model_Theta_tx(Theta),4), sep="", collapse=", " )
  print(msg)

  msg = paste( names(SCNA_model[["mix.w"]]), " = ", round(SCNA_model[["mix.w"]],4), sep="", collapse=", " )
  print(msg)
}

