cutq <- function(x,prefix=character(0)) {
  stopifnot(all(x > 0))
  lbl <- paste0(prefix,paste0("Q", c(1,2,3,4)))
  # get quantiles for non-placebos
  q <- quantile(x,c(0,0.25,0.5,0.75,1)) %>% unname
  # cut non-placebo records
  qs <- cut(x,breaks=q,include.lowest=TRUE,right=FALSE,labels=lbl)
  qs <- as.character(qs)
  # output
  y <- vector(mode="character", length=length(x))
  #' Assign placebo records
  # assign non-placebo records
  y <- qs
  # Order
  y <- factor(y,levels=c(lbl))
  y
}


datag <- function(amt) {
  days <- seq(0,20)
  days <- c(days,28+days)
  out <- vector(mode="list",length=length(amt))
  for(i in seq_along(amt)) {
    out[[i]] <- data_frame(amt=amt[i],time=days,evid=1,cmt=12,ID=i)
  }
  bind_rows(out)
}


.colSet1 <- function(...) ggplot2::scale_color_brewer(palette="Set1",...)
.fillSet1 <- function(...) ggplot2::scale_fill_brewer(palette="Set1",...)

nfact <- function(x,prefix="", suffix="",pad=TRUE) {
  ux <- sort(unique(x))
  if(pad) return(factor(x,ux, paste(prefix,ux,suffix)))
  return(factor(x,ux, paste0(prefix,ux,suffix)))
}
