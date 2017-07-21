
library(mrgsolve)
library(readr)
library(purrrlyr)
library(purrr)
library(dplyr)
library(dmutate)
library(ggplot2)
library(magrittr)
source("functions.R")
noline <- element_blank()
theme_plain <- function(...) {
  theme_bw() + theme(panel.grid.major=noline,panel.grid.minor=noline, 
                     plot.margin=margin(0.5,0.5,1,0.5,unit="cm"),...)
}


mod <- mread("mapk", "model")

source("objects.R")
rotx <- function(angle=30) theme(axis.text.x = element_text(angle = angle, hjust = 1))
roty <- function(angle=30) theme(axis.text.y = element_text(angle = angle, hjust = 1))

data.file <- file.path("data", "s10vpop.csv")

set.seed(1001100)
vp <- read_csv(data.file) %>% sample_n(250,weight=PW,replace=TRUE)

vp %<>% mutate(VPOP2 = 1:n())

cov2 <- cl2+ka2+V2~rlmvnorm(log(c(cl2,ka2,V2)),iiv_vemu) | VPOP2
cov3 <- cl3+ka3+q2+V3+V3b~rlmvnorm(log(c(cl3,ka3,q2,V3,V3b)),iiv_cobi) | VPOP2
cov4 <- cl4+ka4+V4~rlmvnorm(log(c(cl4,ka4,V4)),iiv_gdc) | VPOP2
cov_set <- covset(cov2,cov3,cov4)
vp %<>% mutate_random(cov_set,envir=pki)
vp %<>% mutate(ke2 = cl2/V2,ke3 = cl3/V3,ke4 = cl4/V4)

l <- split(vp,vp$VPOP2)

sim <- function(x,Data) {
  mod %>%
    param(x) %>%
    init(x) %>% Req(TUMOR) %>%
    mrgsim(data=Data,end=56,delta=0.25) %>%
    filter(time==56)
}
siml <- function(data) {
  parallel::mclapply(l, mc.cores=4,function(x) {
    sim(x,Data=data) %>% mutate(VPOP2= x$VPOP2,VPOP=x$VPOP)
  }) %>% bind_rows
}
comb <- function(...) {
  bind_rows(list(...)) %>% arrange(time)
}

# CETUX_dose = 450;   %mg, weekly
# VEMU_dose = 960;    %mg, BID
# COBI_dose = 60;     %mg, daily
# GD ERKI_dose = 400;    %mg, daily 

##' Nothing
data0 <- data_frame(amt=0,evid=1,cmt=8,ID=1,time=0)

##' BFRAF-i CMT 8  - vemurafanib VEMU
dataV <- data_frame(amt=960,evid=1,cmt=8,ID=1,ii=0.5,addl=120,time=0)

##' ERKi CMT 12 - GDC-0994 
dataG <- expand.ev(amt=400,cmt=12,time=c(0,28),ii=1,addl=20) %>% mutate(ID=1)

out <- mrgsim(mod, data=dataG, end=56)
plot(out, ERKi~time)

## MEKi CMT 10 - cobimetinib COBI
dataCO <- mutate(dataG,amt=60,cmt=10)


##' EGFR CMT 7 - cetuximab CETUX
dataCE <- data_frame(time=0,cmt=7,ii=7,addl=7,evid=1,ID=1,amt=450)

comb(dataCE,dataV)


sim1 <- siml(data0 ) %>% mutate(label=1)
sim2 <- siml(dataCE) %>% mutate(label=2)
sim3 <- siml(dataV ) %>% mutate(label=3)
sim4 <- siml(dataCO) %>% mutate(label=4)
sim5 <- siml(dataG ) %>% mutate(label=5)

sim23 <- comb(dataCE, dataV)  %>% siml %>% mutate(label=23)
sim24 <- comb(dataCE, dataCO) %>% siml %>% mutate(label=24)
sim25 <- comb(dataCE, dataG)  %>% siml %>% mutate(label=25)
sim34 <- comb(dataV,  dataCO) %>% siml %>% mutate(label=34)
sim35 <- comb(dataV,  dataG)  %>% siml %>% mutate(label=35)
sim45 <- comb(dataCO, dataG)  %>% siml %>% mutate(label=45)

sim234 <- comb(dataCE, dataV,  dataCO) %>% siml %>% mutate(label=234)
sim235 <- comb(dataCE, dataV,  dataG)  %>% siml %>% mutate(label=235)
sim245 <- comb(dataCE, dataCO, dataG)  %>% siml %>% mutate(label=245)
sim345 <- comb(dataV,  dataCO, dataG)  %>% siml %>% mutate(label=345)


sim2345 <- comb(dataCE,dataV,dataCO,dataG) %>% siml %>% mutate(label=2345)

lab <- c("No TREAT", "CETUX", "VEMU", "COBI", "GDC",
         "CETUX+VEMU", "CETUX+COBI", "CETUX+GDC", "VEMU+COBI","VEMU+GDC", "COBI+GDC",
         "CETUX+VEMU+COBI", "CETUX+VEMU+GDC", "CETUX+COBI+GDC", "VEMU+COBI+GDC",
         "CETUX+VEMU+COBI+GDC")

sims <- bind_rows(sim1,sim2,sim3,sim4,sim5,sim23,sim24,sim25,sim34,
                  sim35,sim45,sim234,sim235,sim245,sim345,sim2345)

ulab <- unique(sims$label)
sims %<>% mutate(labelf = factor(label,levels=ulab,labels=as.character(ulab)))
sims %<>% mutate(labelff = factor(label,levels=ulab,labels=lab))


p1 <- 
  ggplot(data=sims) + 
  geom_point(aes(x=labelff, y=TUMOR),position=position_jitter(width=0.15),col="grey") +
  scale_y_continuous(limits=c(0,2.5),name="Tumor size",breaks=c(0,0.5,1,1.5,2,2.5,3)) +
  scale_x_discrete(name="") + 
  geom_hline(yintercept=0.7,col="firebrick", lty=1,lwd=1)  +
  geom_boxplot(aes(x=labelff, y=TUMOR),fill="darkslateblue",col="darkslateblue",alpha=0.2) +
  theme_plain() + rotx(30)

p1

pdf(file="all_treat_.pdf", width=8,height=4)
p1
dev.off()
