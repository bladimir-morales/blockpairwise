Fisher <- function(par, data){
  n <- nrow(data)
  mu <- rep(par[1],n)
  sigma2 <- par[2]
  phi <- par[3]
  tau2 <- par[4]
  D <-  as.matrix(stats::dist(data[,1:2], method = 'euclidean'))
  R <- exp(-D/phi)
  Sigma <- tau2*diag(n) + sigma2*R

  L <- try(chol(Sigma), silent = TRUE)#Cholesky
  if (inherits(L, "try-error")) {
    Sigma <- Sigma + diag(1e-6, n)
    L <- chol(Sigma)
  }

  solveSigma_mat <- function(L, mat){
    z <- backsolve(L, mat, transpose=T)
    x <- backsolve(L, z)
    return(x)
  }

  I <- matrix(NA,length(par),length(par))
  v1 <- rep(1,n)

  #backsolve(L,(backsolve(L,matrix(v1))))

  dR <- R*(as.matrix(D)/phi^2)

  I[1,1] <- t(v1)%*%solveSigma_mat(L, matrix(v1))
  I[1,2:4] <- 0
  I[2:4,1] <- 0
  SR <- solveSigma_mat(L, R)
  SdR <- solveSigma_mat(L, dR)
  SI <- solveSigma_mat(L, diag(n))

  I[2,2] <- 0.5*sum(diag(SR%*%SR))
  I[2,3] <- 0.5*sigma2*sum(diag(SR%*%SdR))
  I[2,4] <- 0.5*sum(diag(SR%*%SI))
  I[3,2] <- t(I[2,3])
  I[4,2] <- t(I[2,4])

  I[3,3] <- 0.5*sigma2^2*sum(diag(SdR%*%SdR))
  I[3,4] <- 0.5*sigma2*sum(diag(SdR%*%SI))
  I[4,3] <- t(I[3,4])
  I[4,4] <- 0.5*sum(diag(SI%*%SI))
  I <- (I + t(I))/2
  return(I)
}
