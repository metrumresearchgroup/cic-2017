##' ---
##' output: 
##'   md_document:
##'     variant: markdown_github
##' ---

#+ echo=FALSE, message=FALSE
knitr::opts_chunk$set(message=FALSE,fig.path="img/mapk-R-",comment=".",fig.align="center")

library(mrgsolve)
library(readr)
library(magrittr)
library(mrgsolvetk)
library(dplyr)
library(ggplot2)

source("src/functions.R")

mod <- mread("mapk", "model")

vp <- read_csv("data/s10vpop.csv") %>% slice(1)


##' Update parameters and initial conditions
mod %<>% param(vp) %>% init(vp) %>% update(end=56,delta=0.1)

##' ## Create a data set
dataG <- datag(400)

out <- mrgsim(mod,data=dataG,obsonly=TRUE,Req="ERKi_blood,CELLS")
out

plot(out)


##' ## Dose/response
dataG2 <- datag(amt=c(150,200,300,400))
out <- mrgsim(mod,data=dataG2,obsonly=TRUE,Req="ERKi_blood,CELLS")
out

plot(out)



##' ## Sensitivity analysis - `wOR`
.mod <- update(mod,events=as.ev(dataG,keep_id=FALSE))

out <- sens_unif(.mod, n=200, lower = 0.9, upper=1, pars="wOR",Req="ERKi_blood,TUMOR")
out %<>% mutate(wORq = cutq(wOR))

ggplot(out, aes(time,TUMOR,col=wORq,group=ID)) + 
  geom_line() + geom_hline(yintercept=1, lty=2) +
  .colSet1() 


##' ## Sensitivity analysis - `taui4`
##' 
##' - Adding 30% variability to IC50
##' 
.mod <- update(mod,events=as.ev(dataG,keep_id=FALSE))

out <- sens_norm(.mod, n=200, cv=30, pars="taui4",Req="ERKi_blood,TUMOR")

out %<>% mutate(taui4q = cutq(taui4))


ggplot(out, aes(time,TUMOR,col=taui4q,group=ID)) + 
  geom_line() + geom_hline(yintercept=1, lty=2) +
  .colSet1() 


##' ## Explore doses in the `vpop`
set.seed(111)
vp <- read_csv("data/s10vpop.csv") %>% sample_n(250,replace=TRUE,weight=PW)
vp %<>% mutate(ID = 1:n())
data <- datag(1)

sim <- function(dose) {
  data %<>% mutate(amt=dose)
  mrgsim(mod,idata=vp,events=as.ev(data,keep_id=FALSE),
         end=-1,add=56,obsonly=TRUE,Req="TUMOR") %>% mutate(dose=dose)
}

library(parallel)
out <- mclapply(seq(0,400,100),sim)

sims <- bind_rows(out)

ggplot(data=sims, aes(x=factor(dose),y=TUMOR)) + 
  geom_point(position=position_jitter(width=0.1),col="darkgrey") +
  geom_hline(yintercept=0.7,col="firebrick") + 
  geom_boxplot(fill=NA) + ylim(0,2.25)

##' ### Summary
sims %>% 
  group_by(dose) %>%
  summarise(med=median(TUMOR), R30 = mean(TUMOR < 0.7))























