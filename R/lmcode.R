#' Linear Mixed Codes - an extension of Linear Mixed Models
#'
#'
#' @param data_list List of data matrices of varying dimensionality. Attempts to find similarities among all datasets with a core structure.
#' @param config Configuration parameters (required, default provided)
#'
#' @return Main
#'
#' @export
lmcode <- function(data_list,
                   config = lmcode::extract_config(verbose = F)
){

  runtime.start <- Sys.time()

  set.seed(config$seed)

  convergence.parameters <- list(count=0,score.vec=c())

  initialise.model <- initialise.lmcode(data_list = data_list,
                                       config = config)

  main.parameters <- initialise.model$main.parameters
  main.code <- initialise.model$main.code

  if (config$verbose){
    print(paste("Beginning lmcode learning with:      Latent invariant dimension (config$j_dim): ", config$j_dim, "    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }


  while (T){

    prev.code <- main.code

    return_update <- update_set(
      data_list = data_list,
      main.parameters = main.parameters,
      main.code = main.code
    )

    main.parameters <- return_update$main.parameters
    main.code <- return_update$main.code

    total.mse <- mean(abs(main.code$lm.code - prev.code$lm.code))

    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, total.mse)
    MSE <- tail(convergence.parameters$score.vec,1)
    prev.MSE <- tail(convergence.parameters$score.vec,2)[1]

    if (convergence.parameters$count>=1){
      print(paste("Iteration: ",convergence.parameters$count," with Tolerance of: ", abs(prev.MSE - MSE),sep=""))
      if ((convergence.parameters$count >= config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
        break
      }
    }

    convergence.parameters$count = convergence.parameters$count + 1

  }



  runtime.end <- Sys.time()

  return_list <- list(

    main.parameters = main.parameters,

    main.code = main.code,

    meta.parameters = list(
      config = config
    ),

    convergence.parameters = convergence.parameters,

    runtime = list(runtime.start = runtime.start,
                   runtime.end = runtime.end,
                   runtime.total = runtime.end - runtime.start)

  )



  if (config$verbose){
    print(paste("Done! Total runtime of   ", runtime.end - runtime.start ,sep=""))
  }

  return(return_list)

}

#' @export
update_set <- function(data_list,
                       main.parameters,
                       main.code){

  main.parameters$A <- (MASS::ginv(t(main.code$lm.code)%*%(main.code$lm.code))%*%t(main.code$lm.code)%*%data_list$y)
  main.code$lm.code = data_list$y%*%t(main.parameters$A)%*%MASS::ginv(main.parameters$A%*%t(main.parameters$A))

  main.parameters$B <- (MASS::ginv(t(main.code$lm.code)%*%(main.code$lm.code))%*%t(main.code$lm.code)%*%data_list$x)
  main.code$lm.code = data_list$x%*%t(main.parameters$B)%*%MASS::ginv(main.parameters$B%*%t(main.parameters$B))

  main.parameters$beta <- MASS::ginv(t(main.parameters$B)%*%main.parameters$B)%*%t(main.parameters$B)%*%main.parameters$A
  main.parameters$var_y <- (main.parameters$A - main.parameters$B%*%main.parameters$beta)%*%t(main.parameters$A - main.parameters$B%*%main.parameters$beta)

  main.parameters$C <- (MASS::ginv(t(main.code$lm.code)%*%(main.code$lm.code))%*%t(main.code$lm.code)%*%data_list$z)
  main.code$lm.code = data_list$z%*%t(main.parameters$C)%*%MASS::ginv(main.parameters$C%*%t(main.parameters$C))

  main.parameters$D <- diag(diag(main.parameters$var_y - (main.parameters$C%*%main.parameters$u)%*%t(main.parameters$C%*%main.parameters$u)))
  main.parameters$u <- MASS::ginv(t(main.parameters$C)%*%main.parameters$C)%*%t(main.parameters$C)%*%(main.parameters$var_y - main.parameters$D)%*%main.parameters$C%*%MASS::ginv(t(main.parameters$C)%*%main.parameters$C)%*%main.parameters$u%*%MASS::ginv(t(main.parameters$u)%*%main.parameters$u)

  return(list(main.parameters = main.parameters,
              main.code = main.code
  ))

}


pinv <- function(X){
  MASS::ginv(t(X)%*%X)
}


chunk <- function(x,n){
  if (n==1){
    list(x)
  }
  else{
    split(x, cut(seq_along(x), n, labels = FALSE))
  }
}

#' @export
soft_threshold <- function(param,config){

  alpha <- config$regularise$a
  lambda <- config$regularise$l

  gamma <- lambda*alpha

  param <- (((param - gamma)*(param>0)+
               (param + gamma)*(param<0))*(abs(param)>gamma))

  return(param / (1 + lambda*(1 - alpha)))

}

