#' Linear Mixed Forms - an extension of Linear Mixed Models with encoding transformations
#'
#'
#' @param data_list List of data matrices of varying dimensionality. Attempts to find similarities among all datasets with a core structure.
#' @param config Configuration parameters (required, default provided)
#'
#' @return Main
#'
#' @export
lmform <- function(data_list,
                   config = lmform::extract_config(verbose = F)
){

  runtime.start <- Sys.time()

  set.seed(config$seed)

  convergence.parameters <- list(count=0,score.vec=c(),curr.cost = 1e6)

  initialise.model <- initialise.lmform(data_list = data_list,
                                        config = config)

  main.parameters <- initialise.model$main.parameters
  main.form <- initialise.model$main.form

  if (config$verbose){
    print(paste("Beginning lmform learning with:      Latent invariant dimension (config$j_dim): ", config$j_dim, "    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }


  while (T){

    prev.cost <- main.form$K+main.form$H

    for (draw.i in c(1:3)){

      config$type_draw <- draw.i

      return_update <- update_set(
        data_list = data_list,
        main.parameters = main.parameters,
        main.form = main.form,
        config = config
      )

      main.parameters <- return_update$main.parameters
      main.form <- return_update$main.form

    }


    convergence.parameters$curr.cost <- main.form$K+main.form$H
    total.mse <- mean(abs(prev.cost - convergence.parameters$curr.cost))

    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, total.mse)
    MSE <- tail(convergence.parameters$score.vec,1)
    prev.MSE <- tail(convergence.parameters$score.vec,2)[1]

    if (convergence.parameters$count>=1){
      if (config$verbose){
        print(paste("Iteration: ",convergence.parameters$count," with Tolerance of: ", abs(prev.MSE - MSE),sep=""))
      }
      if ((convergence.parameters$count >= config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
        break
      }
    }

    convergence.parameters$count = convergence.parameters$count + 1

  }



  runtime.end <- Sys.time()

  return_list <- list(

    main.parameters = main.parameters,

    main.form = main.form,

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
                       main.form,
                       config){

  main.parameters$A <- (MASS::ginv(t(main.form$K+main.form$H)%*%(main.form$K+main.form$H))%*%t(main.form$K+main.form$H)%*%(data_list$y))
  main.parameters$alpha <- soft_threshold(t(main.parameters$A), config = config)

  main.parameters$B <- (MASS::ginv(t(main.form$K)%*%(main.form$K))%*%t(main.form$K)%*%data_list$x)
  main.parameters$beta <- soft_threshold(t(main.parameters$B), config = config)

  main.parameters$C <- (MASS::ginv(t(main.form$H)%*%(main.form$H))%*%t(main.form$H)%*%data_list$z)
  main.parameters$u <- soft_threshold(t(main.parameters$C), config = config)

  main.form$K <- data_list$x%*%main.parameters$beta%*%MASS::ginv(t(main.parameters$beta)%*%main.parameters$beta)
  main.form$H <- (data_list$y%*%main.parameters$alpha - data_list$x%*%main.parameters$beta)%*%MASS::ginv(t(main.parameters$u)%*%main.parameters$u)

  main.parameters$intercept <- data_list$y%*%main.parameters$alpha - data_list$x%*%main.parameters$beta

  return(list(main.parameters = main.parameters,
              main.form = main.form
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
