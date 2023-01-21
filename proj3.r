#Qunoot Khaleeq - s2314595
# Smoothing and function estimation by combining basis expansions and penalized regression

pspline = function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  # Function to fit P-splines to x, y data, 
  # and choos the smoothing parameter by GCV
  
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  qrx = qr(X)
  R <- qr.R(qrx)
  penalty_matrix <- crossprod(D)
  n <- nrow(X)
  p <- ncol(X)
  I <- diag(p)
  R_inv <- backsolve(R, I) # R^-1
  eigen_decompose = eigen(t(R_inv) %*% (penalty_matrix %*% R_inv))
  eigen_values <- matrix(diag(eigen_decompose$values), ncol = p)
  U <- eigen_decompose$vectors
  R_inv_U <- R_inv %*% U
  UT_QTY <- crossprod(U, (qr.qty(qrx, y)[1:p]))
  log_lambdas <- seq(logsp[1], logsp[2], length = ngrid)
  gcv <- c()
  for(i in log_lambdas){
    inversed <- solve(I + (exp(i)*eigen_values))
    beta_hat <- R_inv_U %*% (inversed %*% UT_QTY)
    edf = sum(diag(inversed))
    sigma_square <- sum((y - (X %*% beta_hat))^2) / (n-edf)
    gcv <- c(gcv, sigma_square/(n-edf))
  }
  
  lambda = exp(log_lambdas[max(which(gcv == min(gcv)))])
  inversed = solve(I + (lambda*eigen_values))
  prod <- R_inv_U %*% inversed
  beta_hat <- prod %*% UT_QTY
  edf = sum(diag(inversed))
  Mu = X %*% beta_hat
  sigma_square <- sum((y - Mu)^2) / (n-edf)
  beta <- (prod %*% crossprod(U, t(R_inv))) * sigma_square # covariance
  
  sumarise <- list(coef=beta_hat, fitted=Mu, sig2=sigma_square, 
              knots=knots, bord=bord, pord=pord, edf=edf,
              x=x, y=y, n=n, V=beta)
  
  class(sumarise) <- "pspline"
  sumarise
}

print.pspline <- function(m){
  
  # print the summary of object of class pspline
  
  cat(sprintf("Order %d p-spline with order %d penalty\n", 
              m$bord, m$pord))
  cat(sprintf("Effective degrees of freedom: %f\b\tCoefficients: %d\n", 
              m$edf, nrow(m$coef)))
 
  r2 <- 1 - ((m$n − 1) * m$sig2)/sum((m$y-mean(m$y))^2)
  gcv <- m$sig2/(m$n - m$edf)
  
  cat(sprintf("residual std dev: %f\tr-squared: %f\tGCV: %f\n", 
              sqrt(m$sig2), r2, gcv))
  
  elements = list(gcv = gcv, edf = m$edf, r2 = r2)
}

predict.pspline <- function(m,x,se=TRUE) {
  Xp <- splines::splineDesign(m$knots,x,ord=m$bord+1,outer.ok=TRUE)
  fit <- Xp %*% m$coef
  if(se){
    std_error <- rowSums(Xp * (Xp %*% m$V))^0.5
    predictions <- list(fit=fit, se=std_error)
  }
  else{
    predictions <- fit
  }
  predictions
}

plot.pspline  <- function(m){
  # Plot graphs
  plot(m$x, m$y, main="Raw Data")
  plot((m$y-m$fitted), m$fitted, main="Residuals vs fitted values")
  qqnorm(m$y-m$fitted, pch = 1, frame = FALSE)
} 
