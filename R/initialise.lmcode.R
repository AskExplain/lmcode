#' @export
initialise.lmcode <- function(data_list,
                             config
                             ){


  main.code <- list(lm.code = array(rnorm(dim(data_list$x)[1]*config$j_dim),dim=c(dim(data_list$x)[1],config$j_dim)))

  # main.parameters <- list(alpha = irlba::irlba(data_list$y,nv = config$j_dim,maxit = 5)$v,
  #                         beta = irlba::irlba(data_list$x,nv = config$j_dim,maxit = 5)$v,
  #                         u = irlba::irlba(data_list$z,nv = config$j_dim,maxit = 5)$v
  #                         )

  main.parameters <- list(alpha = array(rnorm(dim(data_list$y)[2]*config$j_dim),dim=c(dim(data_list$y)[2],config$j_dim)),
                          beta = array(rnorm(dim(data_list$x)[2]*config$j_dim),dim=c(dim(data_list$x)[2],config$j_dim)),
                          u = array(rnorm(dim(data_list$z)[2]*config$j_dim),dim=c(dim(data_list$z)[2],config$j_dim))
  )




  return(list(main.parameters=main.parameters,
              main.code=main.code))

}

#' @export
extract_config <- function(verbose=T){
  config <- list(
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
