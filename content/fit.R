
library(readr)
library(dplyr)
library(mrgsolve)
library(minqa)
library(magrittr)
library(tidyr)
library(ggplot2)
source("functions.R")

typef <- function(x) factor(x, c(1,2), c("Pitavastatin alone", "Pitavastatin + CsA"))

data <- read_csv("fig4a.csv") %>% mutate(profile = NULL,type=ID,typef=typef(ID))

ggplot(data=data,aes(time,DV)) + 
  geom_point(col="firebrick") + 
  facet_wrap(~typef) + 
  scale_y_continuous(trans="log", limits=c(0.1,300), breaks=logbr())


dose <- filter(data,evid==1) %>% mutate(typef=NULL)

yobs <- filter(data,evid==0) %>% select(DV) %>% unlist %>% unname

wt <- function(x) 1/x^2


mod <- mread("yoshikado","model") %>% 
  update(end=14,delta=0.1) %>% Req(CP) %>% obsonly

data %<>% select(-typef)

mod %>% 
  mrgsim(data=dose,obsaug=TRUE) %>% 
  plot(CP~.,scales=list(y=list(log=TRUE)))

##' Prediction function
pred <- function(p, .data, yobs=NULL, pred=FALSE) {
  
  .par <- lapply(p,exp) %>% setNames(names(theta))
  .mod <- param(mod,.par)
  
  if(pred) {
    out <- mrgsim(.mod,data=.data,carry.out="type")
    return(as_data_frame(out))
  }
  
  out <- mrgsim(.mod,data=.data,obsonly=TRUE,Req="CP")
  
  return(sum(wt(yobs)*(yobs-out$CP)^2))
  
  #return(-1*sum(dnorm(log(yobs),log(out$CP),.par$sigma,log=TRUE)))
  
}

##' Fit 5 parameters on log scale
theta <- log(c(fbCLintall = 1.2, ikiu = 1.2, 
               fbile = 0.8, ka = 0.1, ktr = 0.1))

control <- list(iprint=25)
fit <- newuoa(theta, pred,.data=data, yobs=yobs,control=control)

fit$par <- setNames(fit$par,names(theta))

df_pred <- pred(fit$par,dose,pred=TRUE) %>% mutate(type = typef(type))
df_init <- pred(theta,dose,pred=TRUE) %>% mutate(type = typef(type))
df_obs <- mutate(data, type=typef(type))


ggplot(data=df_pred) + 
  geom_line(data=df_init,aes(time,CP,lty="A"), col="black", lwd=0.7) +
  geom_line(aes(time,CP,lty="B"),col="darkslateblue",lwd=0.7) + 
  geom_point(data=df_obs,aes(time,DV),col="firebrick",size=2) + 
  facet_wrap(~type) + 
  scale_y_continuous(trans="log",breaks=10^seq(-4,4), 
                     limits=c(0.1,100),
                     "Pitavastatin concentration (ng/mL)") +
  scale_x_continuous(name="Time (hours)", breaks=seq(0,14,2)) +
  scale_linetype_manual(values= c(2,1),labels=c("Initial estimates", "Final estimates"), name="") +
  theme(legend.position="top")

##' OFV value
pred(fit$par,.data=data,yobs=yobs)

exp(fit$par)



library(mrgsolvetk)
library(optimhelp)
ofv <- function(dv,pred,par) {
  -1*sum(dnorm(log(dv),
               log(pred),par$sigma,log=TRUE))
}
ofv <- function(dv,pred,par) {
   sum((1/dv^2)*(dv-pred)^2)
}
par <- parset(log_par("fbCLintall",1.2),
              log_par("ikiu", 1.2), 
              logit_par("fbile", 0.4), 
              log_par("ka", 1),
              log_par("ktr", 0.1))

fitt <- fit_optim(mod,data,pred="CP",ofv=ofv,par=par,method="CG",
                 control=list(trace=10))

coef(fitt$pars) 
fitt$value



##' DEoptim
##' "Performs evolutionary global optimization via the 
##'  Differential Evolution algorithm."
library(RcppDE)

lower <- rep(-6,length(theta)) %>% setNames(names(theta))
upper <- rep(4,length(theta)) %>% setNames(names(theta))

set.seed(330303)
control <- DEoptim.control(NP=10*length(theta), CR=0.925, F=0.85,
                           itermax=100,storepopfrom=0)
fit2 <- DEoptim(fn=pred, lower=lower,upper=upper, control=control,
                .data=data, yobs=yobs)

data.frame(DE = exp(fit2$optim$bestmem),newuoa  = exp(fit$par)) %>% signif(3)


##' DA for the plot
pops <- lapply(fit2$member$storepop, as.data.frame)
hx <- bind_rows(pops)
hx <- mutate(hx, iteration=rep(1:control$itermax,each=control$NP))
hx <- mutate(hx, pop = rep(1:control$NP, time=control$itermax))
hxm <- gather(hx, variable, value, 1:5) %>% mutate(value = exp(value))
best <- as_data_frame(fit2$member$bestmemit) %>% mutate(iteration = 1:control$itermax)
bestm <- gather(best,variable,value,1:5) %>% mutate(value = exp(value))


ggplot(data=hxm) + 
  geom_line(aes(iteration,value,group=pop),col="darkslateblue") + 
  geom_line(data=bestm,aes(iteration,value),col="orange",lwd=1) + 
  scale_y_continuous(trans="log", breaks=10^seq(-4,4), name="Parameter value") + 
  facet_wrap(~variable, ncol=2, scales="free_y") 






