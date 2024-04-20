#!/usr/bin/env Rscript

# get data (read preprocessed copy-number data)
cnv <- readRDS('data-cnv.rds')


# get pathway gene lists
library('qusage')
PID <- read.gmt('data/c2.cp.pid.v7.0.symbols.gmt')

# mutation drivers
mut.dbv3 <- read.table('data/mutation_download_tab.txt', sep ='\t', header = T)
mut.dbv3 <- mut.dbv3[which(mut.dbv3$cancer_type_abbr == 'OV'),]
mut.dbv3.driver_gene <- gsub(" ", "", unlist(strsplit(mut.dbv3$driver_gene, ',')))


# cnv drivers
cnv.dbv3 <- read.table('data/CNV_download_tab.txt', sep ='\t', header = T)
cnv.dbv3 <- cnv.dbv3[which(cnv.dbv3$cancer_type_abbr == 'OV'),]
cnv.dbv3.driver_gene <- gsub(" ", "", unlist(strsplit(cnv.dbv3$driver_gene, ',')))

cnv.dbv3.driver_gene <- c(cnv.dbv3.driver_gene)
common.mut.cnv <- intersect(mut.dbv3.driver_gene, cnv.dbv3.driver_gene)
mut.dbv3.driver_gene <- mut.dbv3.driver_gene[-match(common.mut.cnv, mut.dbv3.driver_gene)]


### read nonlinear model data as output of "model_poi.R" function

poi.nonlin <- readRDS('model-stats-poi-nonlinear.rds')


#### actual thresholding analysis

cnv.means <- rowMeans(cnv, na.rm = T)

set.seed(123L)

reps <- 1000L

aug.names <- function(df) {
  df$name <- rownames(df)
  return (df)
}

results <- NULL
for (j in seq_along(PID)) {
  for (r in seq_len(reps)) {
    
    pw.genes <- PID[[j]][ PID[[j]] %in% rownames(poi.nonlin) ]
    
    cnv.dat <- aug.names( poi.nonlin[ pw.genes[ pw.genes %in% cnv.dbv3.driver_gene ] , ] )
    mut.dat <- aug.names( poi.nonlin[ pw.genes[ pw.genes %in% mut.dbv3.driver_gene ] , ] )
    
    cnv.aug <- rbind( cnv.dat, aug.names(poi.nonlin)[ sample( nrow(poi.nonlin) , length(pw.genes) - nrow(cnv.dat) ) , ] )
    mut.aug <- rbind( mut.dat, aug.names(poi.nonlin)[ sample( nrow(poi.nonlin) , length(pw.genes) - nrow(mut.dat) ) , ] )
    
    results <- rbind( results,
                      data.frame( nu = mean( cnv.means[ cnv.aug$name ] ) ,
                                  rho = sqrt( pmax( 1. - sum( cnv.aug$sum.rsq ) / sum( cnv.aug$sum.rsq0 ) , 0. ) ), from = 'cnv', rep = r ),
                      data.frame( nu = mean( cnv.means[ mut.aug$name ] ),
                                  rho = sqrt( pmax( 1. - sum( mut.aug$sum.rsq ) / sum( mut.aug$sum.rsq0 ) , 0. ) ), from = 'mut', rep = r )
    )
    
  }
}

save(results, file = 'path.thresholding.RData')

get.th <- function(x, y) {
  x.y <- unique(sort(c(x, y)))
  get.f <- function(inds, top) {
    counts <- rep(0L, top)
    for (j in inds)
      counts[j] <- counts[j] + 1L
    return (counts)
  }
  fx <- get.f( match(x, x.y), length(x.y) )
  fy <- get.f( match(y, x.y), length(x.y) )
  obj <- ( sum(fx) - cumsum(fx) ) + ( cumsum(fy) - fy )
  return (x.y[ which.min(obj) ] )
}


x <- results[ results$from == 'mut', ]$nu
y <- results[ results$from == 'cnv', ]$nu
get.th( x, y)




