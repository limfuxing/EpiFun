
# Gamma regression with log-link function
# the highest response variable value is treated as censored
# formula = model specification response ~ independent variables (intercept is included automatically)
# data    = optional data frame where the response and indep vars can be found

censGamma.glm <- function(formula,data=NULL) {
  X     <- model.matrix(formula,data=data)
  y     <- model.frame(formula,data=data)[,1]
  xy    <- as.matrix(cbind(y,X))
  xy    <- na.omit(xy)
  y     <- xy[,1]
  cens  <- y==max(y)
  X     <- as.matrix(xy[,-1])
  fit   <- glm(y~X[,-1],family=Gamma(link='log'))
  p.init <- c(fit$coef,-log(sum(fit$residuals^2)/fit$df.r))

  #nlog-likelihood
  gammacens.nlogl <- function(p,y,X,cens) {

   mu   <- exp(X %*% as.matrix(p[-length(p)]))
   scale<- exp(p[length(p)])
   logl <- (1-cens)*dgamma(y,shape=mu/scale,scale=scale,log=TRUE) + cens*pgamma(y,shape=mu/scale,scale=scale,log=TRUE,lower.tail=FALSE)
   -sum(logl)
  }

  out <- optim(p=p.init,fn=gammacens.nlogl,y=y,cens=cens,X=X,hessian=T)
  table  <- data.frame(est=out$par,SE=sqrt(diag(solve(out$hess))))
  table$z<- table$est/table$SE
  table$pvalue <- 2*pnorm(abs(table$z),lower.tail=FALSE)
  rownames(table) <- c('intercept',colnames(X)[-1],'log(scale) parameter')
  print(paste0('log-likelihood at final solution=',-out$v, ' and ', ifelse(out$conv==0,'converged','not-converged')))
  table
}

# simulate data
x <- rnorm(1000)
z <- rnorm(1000)
y <- rgamma(1000,shape=exp(0.5+0.2*x-0.2*z)/2,scale=2)

# censor anything above 5
y[y>5] <- 5

df<- data.frame(a=x,z=z,y=y)
 

# example: run function and produce output
censGamma.glm(y~a+z, data =df)
