#' @export
initialise.lmcode <- function(data_list,
                             config
                             ){


  main.code <- list(lm.code = array(rnorm(dim(data_list$x)[1]),dim=c(dim(data_list$x)[1],config$j_dim)))

  main.parameters <- list(u = array(rnorm(dim(data_list$z)[2]),dim=c(dim(data_list$z)[2],1)))


  return(list(main.parameters=main.parameters,
              main.code=main.code))

}
