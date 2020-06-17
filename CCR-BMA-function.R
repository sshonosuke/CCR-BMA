#-------------------------------------------------------#
##  Function for improved confidence region 
# Y: (n,2)-matrix of summary statistics (mean)
# S: (n,2)-matrix of summary statistics (variance)
# alpha: significance region 
# grid: the number of grids of boundary of confidence region 

CCR.BMA <- function(Y, S, alpha=0.05, grid=200){
  n <- dim(Y)[1]   # the number of studies
  tr <- function(mat){ sum(diag(mat)) }    # trace function
  X <- c()
  for(i in 1:n){  X <- rbind(X, diag(2))  }
  D <- diag(as.vector(t(S)))
  
  # ordinary least squares estimator 
  OLS <- apply(Y, 2, mean)
  resid <- t(Y)-OLS
  
  # moment estimation of random effects variance-covariance matrix
  mat <- 0
  for(i in 1:n){
    mat <- mat + resid[,i]%*%t(resid[,i]) - diag(S[i,])
  }
  hPsi0 <- mat/n   # naive moment estimator 
  
  # Bias of hPsi0
  BB <- 0
  for(i in 1:n){
    BB <- BB + hPsi0 + diag(S[i,])
  }
  hPsi1 <- hPsi0 + BB/n^2   # bias corecred version
  
  # adjustment
  SVD <- svd(hPsi1)
  Lam <- diag(SVD$d*ifelse(SVD$d>0,1,0))
  hPsi <- t(SVD$u)%*%Lam%*%SVD$u   
  
  # GLS for fixed effects
  yy <- as.vector(t(Y))
  PP <- solve( diag(n)%x%hPsi + D )
  hbeta <- as.vector( solve(t(X)%*%PP%*%X)%*%t(X)%*%PP%*%yy )
  
  ## h.ast
  IH <- 0
  for(i in 1:n){  IH <- IH + solve(hPsi+diag(S[i,]))  }
  V <- solve(IH)
  
  a1 <- c()
  a2 <- array(NA, c(n, n, n))
  a3 <- array(NA, c(n, n, n))
  a4 <- matrix(NA, n, n)
  a5 <- matrix(NA, n, n)
  
  # a1
  for(i in 1:n){
    A <- 0
    for(j in 1:n){
      A <- A + solve(hPsi+diag(S[j,]))%*%(hPsi+diag(S[i,]))%*%solve(hPsi+diag(S[j,]))
    }
    a1[i] <- tr(V%*%A%*%V%*%A)
  }
  
  # a2, a3, a4, a5
  for(i in 1:n){
    for(j in 1:n){
      Bi <- hPsi+diag(S[i,])
      Bj <- hPsi+diag(S[j,])
      a4[i,j] <- tr( V%*%solve(Bi)%*%Bj%*%solve(Bi)%*%Bj%*%solve(Bi) )
      a5[i,j] <- tr( solve(Bi)%*%Bj)*tr(V%*%solve(Bi)%*%Bj%*%solve(Bi) )
      for(k in 1:n){
        Bk <- hPsi+diag(S[k,])
        a2[i,j,k] <- ( tr(V%*%solve(Bi)%*%Bj%*%solve(Bk)) )^2
        a3[i,j,k] <- tr( V%*%solve(Bi)%*%Bj%*%solve(Bk)%*%V%*%solve(Bk)%*%Bj%*%solve(Bi) )
      }
    }
  }
  
  xx <- qchisq(1-alpha, 2)     # quantile of chi-square with 2 degrees of freedom
  val <- (xx/16-3/4)*(sum(a1)+sum(a2))+(xx/16-1/4)*sum(a3)+sum(a4)+sum(a5)
  h.ast <- val/n^2
  
  ## confidence region
  Theta <- seq(0,2*pi, length=grid)
  Lam <- eigen(IH)$values
  HH <- t(eigen(IH)$vectors)
  c <- xx*(1+h.ast)
  tx <- sqrt(c/Lam)*rbind(cos(Theta),sin(Theta))
  Region <- t(hbeta+t(HH)%*%tx)
  
  return(Region)
}




