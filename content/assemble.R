
library(mrgsolve)
library(dmutate)
library(tidyverse)

source("functions.R")

source("objects.R")

mod <- mread("mapk", "model", compile = FALSE)



data.file <- file.path("data", "s10vpop.csv")

set.seed(1001100)
vp <- read_csv(data.file) %>% sample_n(250,weight=PW,replace=TRUE)

vp <-  mutate(vp, VPOP2 = seq(n()))

cov2 <- cl2+ka2+V2~rlmvnorm(log(c(cl2,ka2,V2)),iiv_vemu) | VPOP2
cov3 <- cl3+ka3+q2+V3+V3b~rlmvnorm(log(c(cl3,ka3,q2,V3,V3b)),iiv_cobi) | VPOP2
cov4 <- cl4+ka4+V4~rlmvnorm(log(c(cl4,ka4,V4)),iiv_gdc) | VPOP2
cov_set <- covset(cov2,cov3,cov4)
vp <- mutate_random(vp, cov_set, envir=pki)
vp <-  mutate(vp, ke2 = cl2/V2,ke3 = cl3/V3,ke4 = cl4/V4)

cmts <- intersect(names(vp), names(init(mod)))
as_cmt <- function(x) paste0(x, "_0")
vp <- rename_at(vp, cmts, funs(as_cmt))

select(vp, contains("_0")) %>% summary

saveRDS(file = "data/s10vpop_pk.RDS", vp)



