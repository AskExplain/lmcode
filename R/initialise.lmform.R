#' @export
initialise.lmform <- function(data_list,
                              config
){


  main.form <- list(lm.form = array(rnorm(dim(data_list$y)[1]*config$j_dim),dim=c(dim(data_list$y)[1],config$j_dim)))

  if (config$init == "svd"){
    main.parameters <- list(alpha = irlba::irlba(data_list$y,nv = config$j_dim,maxit = 5)$v,
                            beta = if(!is.null(data_list$x)){irlba::irlba(data_list$x,nv = config$j_dim,maxit = 5)$v}else{NULL},
                            u = if(!is.null(data_list$z)){irlba::irlba(data_list$z,nv = config$j_dim,maxit = 5)$v}else{NULL}
    )
  }

  if (config$init == "rnorm"){
    main.parameters <- list(alpha = array(rnorm(dim(data_list$y)[2]*config$j_dim),dim=c(dim(data_list$y)[2],config$j_dim)),
                            beta = if(!is.null(data_list$x)){array(rnorm(dim(data_list$x)[2]*config$j_dim),dim=c(dim(data_list$x)[2],config$j_dim))}else{NULL},
                            u = if(!is.null(data_list$z)){array(rnorm(dim(data_list$z)[2]*config$j_dim),dim=c(dim(data_list$z)[2],config$j_dim))}else{NULL}
    )
  }

  main.parameters$A <- t(main.parameters$alpha)
  main.parameters$B <- t(main.parameters$beta)
  main.parameters$C <- if(!is.null(data_list$z)){t(main.parameters$u)}else{NULL}

  main.form$lm.form = data_list$y%*%t(main.parameters$A)%*%MASS::ginv(main.parameters$A%*%t(main.parameters$A))
  main.parameters$A <- (MASS::ginv(t(main.form$lm.form)%*%(main.form$lm.form))%*%t(main.form$lm.form)%*%data_list$y)
  main.parameters$alpha <- soft_threshold(t(main.parameters$A), config = config)

  main.form$lm.form = data_list$x%*%t(main.parameters$B)%*%MASS::ginv(main.parameters$B%*%t(main.parameters$B))
  main.parameters$B <- (MASS::ginv(t(main.form$lm.form)%*%(main.form$lm.form))%*%t(main.form$lm.form)%*%data_list$x)
  main.parameters$beta <- soft_threshold(param = t(main.parameters$B), config = config)

  if (!is.null(data_list$z)){
    main.form$lm.form = data_list$z%*%t(main.parameters$C)%*%MASS::ginv(main.parameters$C%*%t(main.parameters$C))
    main.parameters$C <- (MASS::ginv(t(main.form$lm.form)%*%(main.form$lm.form))%*%t(main.form$lm.form)%*%data_list$z)
    main.parameters$u <- soft_threshold(t(main.parameters$C), config = config)
  }

  return(list(main.parameters=main.parameters,
              main.form=main.form))

}

#' @export
extract_config <- function(verbose=T){
  config <- list(
    init = "svd",
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
