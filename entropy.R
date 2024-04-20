##################### genes contribution to pathway scores list and entropy

library('qusage')
PID <- read.gmt('data/c2.cp.pid.v7.0.symbols.gmt')

##--- get gene contribution from contrib model per pathway

contrib <- readRDS('pathway-marker-contrib.rds')
path.nonlin <- readRDS('model-stats-poi-nonlinear-pathway.rds')


normalized <- NULL
compl <- NULL
top <- NULL
prob <- NULL
for (i in 1:length(PID)){
  normalized[[names(PID)[i]]] <- abs(unlist(contrib[names(PID)[i]]) -unlist(path.nonlin[names(PID)[i]]))
  
  p <- unlist(normalized[names(PID)[i]]) / sum(unlist(normalized[names(PID)[i]]))
  z <- p*log(p)
  z[!(p > 0.)] <- 0.
  h <- sum(-z)
  c <- exp(h)
  compl [[names(PID)[i]]] <- c
  prob [[names(PID)[i]]] <- sort(tail(sort(p),5), decreasing =T)
  top[[names(PID)[i]]] <- sort(tail(sort(unlist(normalized[names(PID)[i]])),5), decreasing =T)
  names(top[[names(PID)[i]]]) <- sub(".*\\.", "", names(top[[names(PID)[i]]]))
  names(prob[[names(PID)[i]]]) <- sub(".*\\.", "", names(prob[[names(PID)[i]]]))
}
