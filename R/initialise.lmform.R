#' @export
initialise.lmform <- function(data_list,
                              config
){

  main.form <- list(K=NULL,
                    L=NULL,
                    H=NULL)

  main.fixed <- list(alpha=NULL,
                     beta=NULL)

  main.random <- list(u = NULL,V_inv = NULL, P = NULL, sigma_w = NULL, G = NULL)

  main.parameters <- list(main.fixed=main.fixed,
                          main.random=main.random)


  gcode.config <- gcode::extract_config(F)
  gcode.config$init <- c("rnorm","rnorm")
  gcode.config$i_dim <- config$j_dim
  gcode.config$j_dim <- config$j_dim
  gcode.config$verbose <- F
  join <- gcode::extract_join_framework(F)
  join$alpha <- c(1,2,3)
  join$beta <- c(1,2,3)
  join$code <- c(1,2,3)

  gcode.model <- gcode::gcode(data_list = data_list, config = gcode.config, join = join)

  main.form$K <- t(gcode.model$main.parameters$alpha[[2]])%*%gcode.model$main.code$code[[2]]
  main.form$H <- t(gcode.model$main.parameters$alpha[[3]])%*%gcode.model$main.code$code[[3]]
  main.form$L <- main.form$K+main.form$H

  main.parameters$main.fixed$alpha <- t((MASS::ginv(t(main.form$L)%*%(main.form$L))%*%t(main.form$L)%*%(data_list$y)))
  main.parameters$main.fixed$beta <- t((MASS::ginv(t(main.form$K)%*%(main.form$K))%*%t(main.form$K)%*%(data_list$x)))
  main.parameters$main.random$u <- t(MASS::ginv(t(main.form$H)%*%main.form$H)%*%t(main.form$H)%*%data_list$z)

  main.parameters$main.random$V_inv <- diag(dim(data_list$y)[1]) - main.form$H%*%MASS::ginv(MASS::ginv(t(main.parameters$main.random$u)%*%(main.parameters$main.random$u)) + t(main.form$H)%*%main.form$H)%*%t(main.form$H)
  main.parameters$main.random$P <- main.parameters$main.random$V_inv - main.parameters$main.random$V_inv%*%main.form$H%*%MASS::ginv(t(main.form$H)%*%main.parameters$main.random$V_inv%*%main.form$H)%*%t(main.form$H)%*%main.parameters$main.random$V_inv
  main.parameters$main.random$sigma_w <- (1/dim(data_list$z)[1])*c(t(data_list$y%*%main.parameters$main.fixed$alpha - data_list$x%*%main.parameters$main.fixed$beta - data_list$z%*%main.parameters$main.random$u)%*%(data_list$y%*%main.parameters$main.fixed$alpha - data_list$x%*%main.parameters$main.fixed$beta - data_list$z%*%main.parameters$main.random$u) - sum(diag(main.parameters$main.random$P)))
  main.parameters$main.random$G <- diag(dim(data_list$z)[2])*mean(((1/main.parameters$main.random$sigma_w)*main.parameters$main.random$u%*%t(main.parameters$main.random$u) - t(data_list$z)%*%main.parameters$main.random$P%*%data_list$z))
  main.parameters$main.random$u <- main.parameters$main.random$G%*%t(data_list$z)%*%main.parameters$main.random$V_inv%*%((data_list$y)%*%main.parameters$main.fixed$alpha-(data_list$x)%*%main.parameters$main.fixed$beta)

  main.form$K <- ((data_list$x)%*%main.parameters$main.fixed$beta)%*%MASS::ginv(t(main.parameters$main.fixed$beta)%*%main.parameters$main.fixed$beta)
  main.form$H <- (data_list$z%*%main.parameters$main.random$u)%*%MASS::ginv(t(main.parameters$main.random$u)%*%main.parameters$main.random$u)
  main.form$L <- main.form$H + main.form$K

  main.parameters$main.random$residuals <- data_list$y%*%main.parameters$main.fixed$alpha-data_list$x%*%main.parameters$main.fixed$beta
  main.parameters$main.random$var_Y <- diag(diag((main.parameters$main.random$residuals)%*%t(main.parameters$main.random$residuals)/dim(data_list$y)[1]))
  main.parameters$main.random$D <- diag(diag(main.parameters$main.random$var_Y - data_list$z%*%main.parameters$main.random$G%*%t(data_list$z)))

  return(list(
              main.parameters=main.parameters,
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
