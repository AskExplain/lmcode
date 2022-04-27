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

  data_list <- list(y=data_list$y,x=data_list$x,z=data_list$z)
  gcode.model <- gcode::gcode(data_list = data_list, config = gcode.config, join = join)

  main.form$fixed_form <- (t(gcode.model$main.parameters$alpha[[1]])%*%gcode.model$main.code$code[[1]]+t(gcode.model$main.parameters$alpha[[2]])%*%gcode.model$main.code$code[[2]])/2








  main.parameters$main.fixed$alpha <- t((MASS::ginv(t(main.form$fixed_form)%*%(main.form$fixed_form))%*%t(main.form$fixed_form)%*%(data_list$y)))
  main.parameters$main.fixed$beta <- t((MASS::ginv(t(main.form$fixed_form)%*%(main.form$fixed_form))%*%t(main.form$fixed_form)%*%(data_list$x)))
  main.parameters$main.random$u <- t((MASS::ginv(t(main.form$fixed_form)%*%(main.form$fixed_form))%*%t(main.form$fixed_form)%*%(data_list$z)))

  main.parameters$main.fixed$alpha_projection <- t(main.parameters$main.fixed$alpha)%*%MASS::ginv(main.parameters$main.fixed$alpha%*%t(main.parameters$main.fixed$alpha))
  main.parameters$main.random$residuals <- data_list$y-data_list$x%*%main.parameters$main.fixed$beta%*%main.parameters$main.fixed$alpha_projection
  main.form$random_form <- (data_list$z%*%main.parameters$main.random$u)%*%MASS::ginv(t(main.parameters$main.random$u)%*%main.parameters$main.random$u)

  main.parameters$main.random$var_Y <- diag(diag((main.parameters$main.random$residuals)%*%t(main.parameters$main.random$residuals)/dim(data_list$y)[1]))
  main.parameters$main.random$D <- diag(diag(main.parameters$main.random$var_Y - data_list$z%*%t(data_list$z)))
  diag(main.parameters$main.random$D)[diag(main.parameters$main.random$D)<0] <- 1e-8
  main.parameters$main.random$V_inv <- diag(1/diag(main.parameters$main.random$D)) - diag(1/diag(main.parameters$main.random$D))%*%main.form$random_form%*%MASS::ginv(MASS::ginv(t(main.parameters$main.random$u)%*%(main.parameters$main.random$u)) + t(main.form$random_form)%*%diag(1/diag(main.parameters$main.random$D))%*%main.form$random_form)%*%t(main.form$random_form)%*%diag(1/diag(main.parameters$main.random$D))
  main.parameters$main.random$P <- main.parameters$main.random$V_inv - main.parameters$main.random$V_inv%*%main.form$fixed_form%*%MASS::ginv(t(main.form$fixed_form)%*%main.parameters$main.random$V_inv%*%main.form$fixed_form)%*%t(main.form$fixed_form)%*%main.parameters$main.random$V_inv
  main.parameters$main.random$sigma_w <- (1/dim(data_list$z)[1])*c((t(main.parameters$main.random$residuals - data_list$z%*%main.parameters$main.random$u%*%MASS::ginv(t(main.parameters$main.random$u)%*%(main.parameters$main.random$u))%*%t(main.parameters$main.fixed$alpha))%*%(main.parameters$main.random$residuals - data_list$z%*%main.parameters$main.random$u%*%MASS::ginv(t(main.parameters$main.random$u)%*%(main.parameters$main.random$u))%*%t(main.parameters$main.fixed$alpha)) - sum(diag(main.parameters$main.random$P))))
  main.parameters$main.random$sigma_H <- c(t(data_list$y)%*%main.parameters$main.random$P%*%data_list$y)/(dim(data_list$y)[1])
  main.parameters$main.random$G <- diag(diag((((1/main.parameters$main.random$sigma_w)*main.parameters$main.random$u%*%t(main.parameters$main.random$u) - t(data_list$z)%*%main.parameters$main.random$P%*%data_list$z))))
  diag(main.parameters$main.random$G)[diag(main.parameters$main.random$G)<0] <- 1e-8
  main.parameters$main.random$u <- main.parameters$main.random$G%*%t(data_list$z)%*%main.parameters$main.random$V_inv%*%(main.parameters$main.random$residuals)%*%main.parameters$main.fixed$alpha


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
