#' @export
initialise.lmform <- function(data_list,
                              config
){



  main.form <- list(K=(irlba::irlba(Y%*%t(Y)+X%*%t(X)+Z%*%t(Z),nu = config$j_dim, it = 30)$u),H=NULL)
  main.parameters <- list(A=NULL,alpha=NULL,
                          B=NULL,beta=NULL,
                          C=NULL,u=NULL,
                          intercept=NULL)


  main.parameters$A <- (MASS::ginv(t(main.form$K)%*%(main.form$K))%*%t(main.form$K)%*%(data_list$y))
  main.parameters$alpha <- soft_threshold(t(main.parameters$A), config = config)

  main.parameters$B <- (MASS::ginv(t(main.form$K)%*%(main.form$K))%*%t(main.form$K)%*%data_list$x)
  main.parameters$beta <- soft_threshold(t(main.parameters$B), config = config)

  main.form$K <- data_list$x%*%main.parameters$beta%*%MASS::ginv(t(main.parameters$beta)%*%main.parameters$beta)
  main.form$H <- data_list$y%*%main.parameters$alpha%*%MASS::ginv(t(main.parameters$alpha)%*%main.parameters$alpha) - main.form$K

  main.parameters$C <- (MASS::ginv(t(main.form$H)%*%(main.form$H))%*%t(main.form$H)%*%data_list$z)
  main.parameters$u <- soft_threshold(t(main.parameters$C), config = config)

  main.parameters$intercept <- data_list$y%*%main.parameters$alpha - data_list$x%*%main.parameters$beta

  return(list(main.parameters=main.parameters,
              main.form=main.form))

}

#' @export
extract_config <- function(verbose=T){
  config <- list(
    regularise = list(a=0,l=0),
    j_dim = 30,
    max_iter=10,
    seed = 1,
    tol=1,
    verbose = T
  )

  if (verbose == T){
    print(config)
  }

  return(config)
}
