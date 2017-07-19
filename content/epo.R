library(mrgsolve)
library(dplyr)
library(ggplot2)

mod <- mread("epo","model")


mwf <- ev_days(ev(amt=78*100),days="m,w,f",ii=24*7, addl=3)

mwf <- as.ev(mwf)

##' A single simulation at 100IU/kg TIW x 4 weeks
mod %>% 
  ev(mwf) %>%  zero_re %>%
  mrgsim(end=960,digits=6,delta=0.5) %>% 
  plot


##' Once weekly same as TIW
qw <- ev(amt=40000, ii=24*7, addl=3)

data <- as_data_set(mwf,qw)

data

mod %>% 
  data_set(data) %>% zero_re %>%
  mrgsim(end=700,digits=6,delta=0.5) %>%
  plot


##' Tolerance / rebound in reticulocytes
##' What controls the depth?
e2 <- ev_days(ev(amt=70*100),days="m,w,f",ii=24*7, addl=2)

out <- 
  mod %>% 
  ev(mwf) %>% 
  mrgsim(end=2400,digits=6,delta=0.25)
plot(out)

library(mrgsolvetk)


.mod <- update(mod, events=mwf, 
               end=2400,tscale=1/24/7, delta=4) %>% zero_re

##' SENSITIVITY ANALYSIS
##' EPO CLEARANCE
out <- sens_norm(.mod, cv=50,  pars="THETA4",n=50)

out <- tidyr::gather(out,variable,value,c(EPOi:HGBi))

ggplot(out, aes(time,value,group=.n,col=THETA4)) + 
  geom_line() + facet_wrap(~variable, scales="free_y")

##" SENSITIVITY ANALYSIS
##' KD
out <- sens_norm(.mod, cv=50,  pars="THETA17", n=50)

out <- tidyr::gather(out,variable,value,c(EPOi:HGBi))

ggplot(out, aes(time,value,group=.n,col=THETA17)) + 
  geom_line() + facet_wrap(~variable, scales="free_y")


