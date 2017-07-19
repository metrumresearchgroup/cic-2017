iiv <- list()

## RTK
cetux <- list(V1 = 4.5, cl1=0.67)
cetux$ke1 <- with(cetux,cl1/V1)
iiv_1 <- c(0,0,0)
iiv$cetux <- as_bmat(iiv_1)

# RAF / BRAF
vemu <- list(cl2=32.2,ka2=4.48,V2=114)
vemu$ke2 <- with(vemu,cl2/V2)
iiv_2<- c(0.109,0.000, 1.02, 0.106, 0.00, 0.457)
iiv$vemu <- as_bmat(iiv_2)


# MEK
cobi <- list(cl3=327,ka3=33.7,q2=252,V3=487,V3b=335)
cobi$ke3 <- with(cobi,cl3/V3)

iiv_3 <- c(
  0.35,
  0.000, 2.81,
  0.000, 0.00, 0.517,
  0.264, 0.00, 0.000, 0.2600,
  0.180 ,0.00, 0.000, 0.0625, 0.557
)
iiv$cobi <- as_bmat(iiv_3)

# ERK 
gdc <- list(cl4=161,ka4=35,V4=171)
gdc$ke4 <- with(gdc,cl4/V4)
iiv_4 <- c(
  0.182, 
  0.000, 3.000,
  0.000, 0.000, 0.149 
)
iiv$gdc <- as_bmat(iiv_4)

pki <- c(vemu,cobi,gdc,cetux)
pki$iiv_cetux <- iiv$cetux
pki$iiv_vemu <- iiv$vemu
pki$iiv_cobi <- iiv$cobi
pki$iiv_gdc <- iiv$gdc



evcet <- ev(amt=100,ii=2*7,addl=1,cmt=7)
evv <- ev(amt=960,ii=0.5, addl=49,cmt=8)
evcob <- ev(amt=1, ii=1, addl=49,cmt=10)
evg <- ev(amt=10*70,ii=1,addl=1000,cmt=12)


