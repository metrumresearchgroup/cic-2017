##' ---
##' title: ""
##' output: pdf_document
##' output_file: "run_stat.pdf"
##' ---

##' ## Quantitative Analyses of Hepatic OATP-Mediated Interactions 
##' Between Statins and Inhibitors Using PBPK Modeling With a Parameter 
##' Optimization Method

#+ echo=FALSE, message=FALSE
library(mrgsolve)
library(dplyr)
library(ggplot2)
options(mrgsolve_mread_quiet=TRUE)




##' stat model
mod <- mread("yoshikado", "model")

##' CsA simulation
ec <- ev(amt=2000,cmt=2)
out <- mod %>% ev(ec) %>% mrgsim(delta=0.01,end=12)
sims <- as.tbl(out)
ggplot(sims, aes(x=time,y=CSA)) + 
  scale_y_continuous(trans="log", 
                     breaks=c(50,100,500,1000,10000,5000,50000), 
                     limit=c(50,50000)) + 
  scale_x_continuous(breaks=seq(0,12,2)) + 
  geom_line() + 
  geom_line(aes(x=time,y=CSAliv),lty=2) + 
  geom_hline(yintercept=2000)


##' # Pitavastatin (EHC model) simulation
ep <- ev(amt=30, time=1)
out <- mod %>% ev(ep) %>% mrgsim(delta=0.1,add=seq(1,1.5,0.01), end=14)

simsp <- as.tbl(out)
ggplot(simsp, aes(x=time,y=CP)) + 
  scale_y_continuous(trans="log10", limit=c(0.1,100),breaks=c(0.1,1,10,100)) + 
  geom_line() + scale_x_continuous(breaks=seq(0,14,2),limit=c(0,14)) 


out <- mod %>% ev(ep+ec) %>% mrgsim(delta=0.01,end=14)

##' 
##' \newpage
##' 
##' # Blue: pitavastatin, no ddi; red: pitavastatin with CsA ddi
simspc <- as.tbl(out)
ggplot(simspc, aes(x=time,y=CP)) + 
  scale_y_continuous(trans="log10", limit=c(0.1,100), breaks=c(0.1,1,10,100)) + 
  geom_line(col="firebrick") + scale_x_continuous(breaks=seq(0,14,2),limit=c(0,14)) + 
  geom_line(data=simsp, aes(x=time,y=CP),col="darkslateblue")


