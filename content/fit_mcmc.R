
library(readr)
library(dplyr)
library(mrgsolve)
library(minqa)

library(tidyr)
library(ggplot2)

typef <- function(x) factor(x, c(1,2), c("Pitavastatin alone", "Pitavastatin + CsA"))

data.file <- file.path("data", "fig4a.csv")

data <- read_csv(data.file) %>% mutate(profile = NULL,type=ID)

dose <- filter(data,evid==1)

yobs <- filter(data,evid==0) %>% dplyr::select(DV) %>% unlist %>% unname

mod <- mread("yoshikado", "model") %>% 
  update(end=14,delta=0.1) %>% Req(CP) %>% obsonly


##' Log prior density without constants:
nprior <- function(theta,mu=0,tau2=1E-6) {
  -0.5*tau2*(theta-mu)^2
}
igprior <- function(theta,a=0.001,b=0.001) {
  -(a+1)*log(theta) - b/theta
}


##' Returns log prior + log likelihood
mcfun <- function(par,.data,n,yobs,pred=FALSE) {
  
  par <- setNames(par,n)
  
  mod <-  param(mod,lapply(par[which_pk],exp))
  
  out <- mrgsim(mod,obsonly=TRUE,Req="CP",data=.data) 
  
  if(pred) return(out)
  
  if(any(out$CP <= 0)) return(Inf)
  
  log.yhat <- log(out$CP)
  log.y    <- log(yobs)
  
  sig2 <- exp(par[which_sig])
  
  data.like <- dnorm(log.y, 
                     mean = log.yhat, 
                     sd   = sqrt(sig2), 
                     log  = TRUE)
  
  pri.pkpars <- nprior(par[which_pk])
  pri.sig2 <- igprior(sig2)
  jac.sig2 <- log(sig2)
  sum.prior <- sum(pri.pkpars,pri.sig2,jac.sig2)
  
  return(sum(data.like,sum.prior))
}

##' Fit 4 parameters on log scale
theta <- log(c(fbCLintall = 1.2, ikiu = 0.1, 
               fbile = 0.8,  ktr =1, sigma=1))

which_pk <- which(names(theta) != "sigma")
which_sig <- which(names(theta) == "sigma")

library(MCMCpack)
contr <- list(fnscale = -1, trace = 1)
set.seed(1110)
fit <- MCMCmetrop1R(mcfun,theta.init=theta, .data=data, 
                    optim.method="BFGS", verbose=100,
                    n=names(theta), yobs=yobs, tune = 1,
                    optim.control=contr, burnin=5000, mcmc=5000)

plot(fit)

summary(exp(fit))

