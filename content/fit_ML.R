
library(readr)
library(dplyr)
library(mrgsolve)
library(minqa)

library(tidyr)
library(ggplot2)

typef <- function(x) factor(x, c(1,2), c("Pitavastatin alone", "Pitavastatin + CsA"))

data <- read_csv("fig4a.csv") %>% mutate(profile = NULL,type=ID)

dose <- filter(data,evid==1)

yobs <- filter(data,evid==0) %>% select(DV) %>% unlist %>% unname
wt <- function(x) 1/x^2


mod <- mread("yoshikado", "model") %>% 
  update(end=14,delta=0.1) %>% Req(CP) %>% obsonly


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
  
  #return(sum(wt(yobs)*(yobs-out$CP)^2))
  
  if(any(out$CP <= 0)) return(Inf)
  
  return(-1*sum(dnorm(log(yobs),log(out$CP),.par$sigma,log=TRUE)))
  
}

##' Fit 5 parameters on log scale
theta <- log(c(sigma=1, fbCLintall = 1.2, ikiu = 1.2, 
               fbile = 0.8, ka = 0.1, ktr = 0.1))

fit <- newuoa(theta, pred,.data=data, yobs=yobs)

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
