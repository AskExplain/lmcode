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

  convergence.parameters <- list(count=0,score.vec=c(1e6))

  initialise.model <- initialise.lmform(data_list = data_list,
                                        config = config)

  main.parameters <- initialise.model$main.parameters
  main.form <- initialise.model$main.form

  if (config$verbose){
    print(paste("Beginning lmform learning with:      Latent invariant dimension (config$j_dim): ", config$j_dim, "    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }


  while (T){
    score.prev <- abs(main.parameters$main.random$residuals)

    return_update <- update_set(
      data_list = data_list,
      main.parameters = main.parameters,
      main.form = main.form,
      config = config
    )

    main.parameters <- return_update$main.parameters
    main.form <- return_update$main.form


    score.curr <- abs(main.parameters$main.random$residuals)

    MSE <- mean(abs(score.curr-score.prev))
    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, MSE)

    s1 <- tail(convergence.parameters$score.vec,1)
    s2 <- tail(convergence.parameters$score.vec,2)[1]

    if (convergence.parameters$count>=1){
      if (config$verbose){
        print(paste("Iteration: ",convergence.parameters$count," with Tolerance of: ", abs(s2-s1),sep=""))
      }
      if ((convergence.parameters$count >= config$max_iter ) | abs(s2-s1) < config$tol){
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

  main.parameters$main.fixed$alpha <- t((MASS::ginv(t(main.parameters$main.random$V_inv%*%main.form$fixed_form)%*%(main.parameters$main.random$V_inv%*%main.form$fixed_form))%*%t(main.parameters$main.random$V_inv%*%main.form$fixed_form)%*%main.parameters$main.random$V_inv%*%(data_list$y)))
  main.parameters$main.fixed$beta <- t((MASS::ginv(t(main.parameters$main.random$V_inv%*%main.form$fixed_form)%*%(main.parameters$main.random$V_inv%*%main.form$fixed_form))%*%t(main.parameters$main.random$V_inv%*%main.form$fixed_form)%*%main.parameters$main.random$V_inv%*%(data_list$x)))
  main.parameters$main.fixed$alpha_projection <- t(main.parameters$main.fixed$alpha)%*%MASS::ginv(main.parameters$main.fixed$alpha%*%t(main.parameters$main.fixed$alpha))
  main.parameters$main.random$residuals <- data_list$y-data_list$x%*%main.parameters$main.fixed$beta%*%main.parameters$main.fixed$alpha_projection

  main.form$random_form <- (data_list$z%*%main.parameters$main.random$u)%*%MASS::ginv(t(main.parameters$main.random$u)%*%main.parameters$main.random$u)
  main.form$fixed_form <- (data_list$y%*%main.parameters$main.fixed$alpha)%*%MASS::ginv(t(main.parameters$main.fixed$alpha)%*%main.parameters$main.fixed$alpha)

  main.parameters$main.random$var_Y <- diag(diag((main.parameters$main.random$residuals)%*%t(main.parameters$main.random$residuals))/dim(data_list$y)[1])

  main.parameters$main.random$D <- diag(diag(main.parameters$main.random$var_Y - data_list$z%*%main.parameters$main.random$G%*%t(data_list$z)))
  diag(main.parameters$main.random$D)[diag(main.parameters$main.random$D)<0] <- 1e-5

  main.parameters$main.random$V_inv <- diag(1/diag(main.parameters$main.random$D)) - diag(1/diag(main.parameters$main.random$D))%*%main.form$random_form%*%MASS::ginv(MASS::ginv(t(main.parameters$main.random$u)%*%main.parameters$main.random$G%*%(main.parameters$main.random$u)) + t(main.form$random_form)%*%diag(1/diag(main.parameters$main.random$D))%*%main.form$random_form)%*%t(main.form$random_form)%*%diag(1/diag(main.parameters$main.random$D))
  main.parameters$main.random$P <- main.parameters$main.random$V_inv - main.parameters$main.random$V_inv%*%data_list$x%*%MASS::ginv(t(data_list$x)%*%main.parameters$main.random$V_inv%*%data_list$x)%*%t(data_list$x)%*%main.parameters$main.random$V_inv

  main.parameters$main.random$sigma_w <- main.parameters$main.random$sigma_w+(1/dim(data_list$z)[1])*c((t(main.parameters$main.random$residuals - data_list$z%*%main.parameters$main.random$u%*%main.parameters$main.fixed$alpha_projection)%*%(main.parameters$main.random$residuals - data_list$z%*%main.parameters$main.random$u%*%MASS::ginv(t(main.parameters$main.random$u)%*%(main.parameters$main.random$u))%*%t(main.parameters$main.fixed$alpha)) - main.parameters$main.random$sigma_w*sum(diag(main.parameters$main.random$P))))
  main.parameters$main.random$sigma_H <- c(t(data_list$y)%*%main.parameters$main.random$P%*%data_list$y)/(dim(data_list$y)[1])

  main.parameters$main.random$G <- main.parameters$main.random$G+(1/dim(data_list$z)[2])*diag(diag((((1/main.parameters$main.random$sigma_w)*(main.parameters$main.random$u)%*%t(main.parameters$main.random$u) - main.parameters$main.random$G%*%t(data_list$z)%*%main.parameters$main.random$P%*%data_list$z%*%main.parameters$main.random$G))))
  diag(main.parameters$main.random$G)[diag(main.parameters$main.random$G)<0] <- 1e-5

  main.parameters$main.random$u <- main.parameters$main.random$G%*%t(data_list$z)%*%main.parameters$main.random$V_inv%*%(main.parameters$main.random$residuals)%*%main.parameters$main.fixed$alpha

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
