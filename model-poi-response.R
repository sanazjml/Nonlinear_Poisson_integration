#!/usr/bin/env Rscript

source('poilin.R', chdir = T)

# data mangling section
# ----


# get data (read preprocessed cop-number and gene expression data)
cnv <- readRDS('data-cnv.rds')
rna <- readRDS('data-rna.rds')

# match genes
rows <- match( rownames(cnv), rownames(rna) )
rows <- data.frame( cnv.row = seq_along(rows), rna.row = rows )[which(!is.na(rows)), , drop = F]


# match samples
cnv.samps <- unique( sub( '[.]([^.]*)$', '', colnames(cnv) ) )
rna.samps <- colnames(rna)
samp.filter <- function(samp) {
  # drop DNAx/RNAy
  samp <- sub( '_(DNA|RNA)\\d*', '', samp )
  # kill cell lines etc.
  samp <- sub( '^(.*)_CL\\d*(.*)$', '', samp )
  samp <- sub( '^(.*)_xO(.*)$', '', samp )
  samp <- sub( '^(.*)_sph(.*)$', '', samp )
  # drop Datax from RNA samples
  samp <- sub( '_(Data)\\d*'   , '' , samp)
  return (samp)
}
common.samps.filtered <- setdiff( intersect(
  samp.filter(cnv.samps), samp.filter(rna.samps) ), '' )
cols <- data.frame(
  cnv.samp = cnv.samps[ match(common.samps.filtered, samp.filter(cnv.samps)) ],
  rna.samp = rna.samps[ match(common.samps.filtered, samp.filter(rna.samps)) ] )
cols$cnv.maj.col <- match( sprintf('%s.%s', cols$cnv.samp, 'nMajor'), colnames(cnv) )
cols$cnv.min.col <- match( sprintf('%s.%s', cols$cnv.samp, 'nMinor'), colnames(cnv) )
cols$rna.col <- match( cols$rna.samp, colnames(rna) )

# just primary samples
cols <- cols[grep('_p', cols$cnv.samp),]


# loop through
outs <- NULL
for (i in seq_len(nrow(rows))) {
  
  # data section
  # ----
  
  # grab data
  expr.gains <- t(attr(rna, 'gain'))[ cols$rna.col, ]
  data <- data.frame( expr = rna[rows$rna.row[i], cols$rna.col] / expr.gains,
                      maj = cnv[rows$cnv.row[i], cols$cnv.maj.col], min = cnv[rows$cnv.row[i], cols$cnv.min.col] )
  data$tot <- data$maj + data$min   # total cnv
  
  #--- read clinical data and extract response parameters of your interest such as PFI or PFS
  root <- ''
  clinical <- read.table(file.path(root, ''), sep = ',', header =T, row.names = 1L, stringsAsFactors = FALSE)
  
  
  pfi.data <- data.frame( is.prog = clinical$Progression.Yes_No == 'Yes',           # cancer progression recorded
                          time.prog = clinical$PFI.after.Primary.therapy.if.Prog_Days  / days.per.mon ,   # time to progression
                          time.last = clinical$PFI.at.outcome.update.when.no.prog_Days / days.per.mon ,   # time to last follow up
                          row.names = clinical$Patient.card..Patient.cohort.code_Patient.Card )  
  
  pfi.data$is.poor <- with( data.frame( pfi.data, short = 6 ),
                            ifelse( pfi.ub <= short, TRUE,     # PFI <= pfi.ub <= short -- short for sure
                                    ifelse( pfi.lb > short, FALSE,   # PFI >= pfi.lb > short -- not short for sure
                                            NA ) ) )                         # not sure
  
  pfi.data$is.good <- with( data.frame( pfi.data, long = 12 ),
                            ifelse( pfi.lb > long, TRUE,
                                    ifelse( pfi.ub <= long, FALSE,
                                            NA ) ) )
  
  pfi.data$group <- ifelse( pfi.data$is.poor, 'poor',
                            ifelse( pfi.data$is.good, 'good', NA ) )
  
  # ----
  
  
  cols$PFI = pfi.data$group[match(sub("_.*", "", cols$cnv.samp), rownames(pfi.data))]
  
  data$response <- NA
  data$response[which(cols$PFI == "poor") ] <- "poor"
  data$response[which(cols$PFI == "good") ] <- "good"
  
  
  # prune missing values
  mask <- !( rowSums(is.na(data)) > 0 )
  expr.gains <- expr.gains[mask]
  data <- data[mask, ]
  
  data <- rbind(data[grep('poor', data$response),], data[grep('good', data$response),])
  # model section
  # ----
  
  #' fit weighted Poisson models
  poilin.fit.wts <- function(design, resp, weights) {
    # fit
    res <- poilin.fit.sing(design * weights, resp * weights)
    
    # TODO: we should use some mutual information analog of 
    # the correlation to be consistent
    
    # compute also sum-of-squares for correlations
    x <- res$x
    x[is.na(x)] <- 0.
    res$sum.rsq <- colSums( weights * ( as.matrix(design) %*% x - resp )^2 )
    
    return (res)
  }
  
  
  # build up a nonlinear model
  
  top <- max( data$maj )
  design <- outer( data$maj, 0:top, `>=` ) + outer( data$min, 0:top, `>=` )
  rownames(design) <- rownames(data)
  colnames(design) <- sprintf('to.%d', 0:top)    # increase from X-1 to X
  colnames(design)[1] <- 'bas'                   # basal effect (cn=0)
  
  
  # nonlinear model with interactions
  fit.0 <- poilin.fit.wts( design, data.frame( expr = data$expr ), weights = expr.gains )
  
  
  ####### first case interval is mirrored in the block matrix
  design <- outer( data$maj, 0:top, `>=` ) + outer( data$min, 0:top, `>=` )
  
  ## poor block
  design.poor <- design
  design.poor[which(data$response != 'poor'),] =  0
  
  ## good block
  design.good <- design
  design.good[which(data$response != 'good'),] = 0
  
  design.pg <- cbind(design.poor, design.good)
  rownames(design.pg) <- rownames(data)
  colnames(design.pg) <- c('bas.p', sprintf('p.to.%d', 1:top), 'bas.g', sprintf('g.to.%d', 1:top)) 
  
  ## mirroring good block
  row.g = grep('good', data$response)
  design.pg [ dplyr::first(row.g): dplyr::last(row.g), 1:(top+1) ] = design.pg[ dplyr::first(row.g): dplyr::last(row.g) , (top+2): (2*(top+1))]
  
  
  ####### second case poor is mirrored in the block matrix
  
  design.gp <- cbind(design.poor, design.good)
  rownames(design.gp) <- rownames(data)
  colnames(design.gp) <- c('bas.p', sprintf('p.to.%d', 1:top), 'bas.g', sprintf('i.to.%d', 1:top)) 
  
  ## mirroring primary block
  row.p = grep('poor', data$response)
  design.gp [ dplyr::first(row.p): dplyr::last(row.p), (top+2): (2*(top+1)) ] = design.gp[ dplyr::first(row.p): dplyr::last(row.p) , 1:(top+1)]
  
  # function to rewrite stacked parameters into unstacked form
  library('MASS')
  rewrite.coefs <- function(x, A0, A1) `names<-`( ginv(A1) %*% A0 %*% x, colnames(A1) ) 
  
  
  
  fit.pg <- poilin.fit.wts( design.pg, data.frame( expr = data$expr ), weights = expr.gains )
  fit.pg$x[is.na(fit.pg$x)] <- 0.
  fit.pg$old.x <- rewrite.coefs( fit.pg$x, design.pg, design )    # rewrite stacked parameters into unstacked form
  names(fit.pg$old.x) <- c('bas.p', sprintf('p.to.%d', 1:top))
  
  fit.gp <- poilin.fit.wts( design.gp, data.frame( expr = data$expr ), weights = expr.gains )
  fit.gp$x[is.na(fit.gp$x)] <- 0.
  fit.gp$old.x <- rewrite.coefs( fit.gp$x, design.gp, design )    # rewrite stacked parameters into unstacked form
  names(fit.gp$old.x) <- c('bas.g', sprintf('g.to.%d', 1:top))
  
  fit <- fit.pg
  fit$sign <- "good"    #good is higher
  if (fit.gp$sum.res < fit.pg$sum.res){
    fit <- fit.gp
    fit$sign <- "poor"   #poor is higher
  }
  if (fit.gp$sum.res == fit.pg$sum.res)  
    fit$sign <- "none"  
  
  
  
  # test any effects
  va.any <- poilin.select( fit.0, fit )
  
  
  # results section
  # ----
  
  # set up result row
  out.row <- data.frame(samps = nrow(data))    # number of data points used in fit
  rownames(out.row) <- rownames(rna)[ rows$rna.row[i] ]    # gene name
  
  # put in constant model
  out.row$coef.avg <- fit.0$x[1]               # avg expr level from model
  out.row$rs.avg <- fit.0$sum.res           # residual (no longer SS, but LL equivalent) after avg
  
  # put in the full model
  out.row <- cbind( out.row,                   # as commented lines above, but works with the nonlinear model as well
                    `colnames<-`( as.data.frame( t(fit$old.x) ), sprintf('%s.%s', 'coef', names(fit$old.x)) ) )
  out.row$rs.fit <- fit$sum.res              # residual SS after model
  out.row$pv.fit <- va.any$p.value              # p-value for model
  out.row$xv.fit <- va.any$stat
  out.row$rho.fit <- sqrt(pmax( 1. - fit$sum.rsq    / fit.0$sum.rsq , 0. ))   # (absolute value of) correlation
  #                             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ fraction of explained variance
  #                                                   ^^^^^^^^^^^^^^ total variance
  out.row$rho.fit.u <- sqrt( 1. - exp( -2.*pmax( out.row$rs.avg - out.row$rs.fit, 0. ) / ( 2.*out.row$rs.avg ) ) )
  out.row$sum.rsq <- fit$sum.rsq
  out.row$sum.rsq0 <- fit.0$sum.rsq
  
  out.row$sign.fit <- fit$sign
  
  
  # append
  outs[[i]] <- out.row
  
}


# concatenate outputs
source('toporbind.R')
out <- do.call( toporbind, outs )

#Normalized X2 value Data
range01 <- function(x, na.rm = T){(x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))}
normalized <- range01(out$xv.fit, na.rm = T)
out$xv.fit.01 <- normalized


# dump binary
saveRDS(out, 'model-stats-poi-nonlinear-response-primary-u.rds')

# dump text
#write.gz.table <- function(data, filename, ...) {
#	file <- gzfile(filename, 'wb')
#	write.table(data, file, ...)
#	close(file)
#}
#write.gz.table(out, 'model-stats.tsv.gz', sep = '\t', row.names = T, col.names = NA)
