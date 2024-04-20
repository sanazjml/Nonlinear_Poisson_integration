#!/usr/bin/env Rscript

source('poilin.R', chdir = T)


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

	# prune missing values
	mask <- !( rowSums(is.na(data)) > 0 )
	expr.gains <- expr.gains[mask]
	data <- data[mask, ]

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

	# constant model (null model)
	fit.0 <- poilin.fit.wts( data.frame( bas = rep(1., nrow(data)) ), data.frame( expr = data$expr ), weights = expr.gains )

	# linear model with no interactions (null for testing asymmetric)
	fit.lin <- poilin.fit.wts( data.frame( bas = 1., tot = data$tot ), data.frame( expr = data$expr ), weights = expr.gains )

	
	# build up a nonlinear model
	top <- max( data$maj )
	design <- outer( data$maj, 0:top, `>=` ) + outer( data$min, 0:top, `>=` )
	rownames(design) <- rownames(data)
	colnames(design) <- sprintf('to.%d', 0:top)    # increase from X-1 to X
	colnames(design)[1] <- 'bas'                   # basal effect (cn=0)

	# nonlinear model with interactions
	fit <- poilin.fit.wts( design, data.frame( expr = data$expr ), weights = expr.gains )

# test section
# ----

	# test any effects
	va.any <- poilin.select( fit.0, fit )
	# test asymmetric effects (major vs minor)
	va.asy <- poilin.select( fit.lin, fit )

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
		`colnames<-`( as.data.frame( t(fit$x) ), sprintf('%s.%s', 'coef', names(fit$x)) ) )
	out.row$rs.fit <- fit$sum.res              # residual SS after model
	out.row$pv.fit <- va.any$p.value              # p-value for model
	out.row$xv.fit <- va.any$stat
	out.row$rho.fit <- sqrt(pmax( 1. - fit$sum.rsq    / fit.0$sum.rsq , 0. ))   # (absolute value of) correlation
	#                             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ fraction of explained variance
	#                                                   ^^^^^^^^^^^^^^ total variance
	out.row$rho.fit.u <- sqrt( 1. - exp( -2.*pmax( out.row$rs.avg - out.row$rs.fit, 0. ) / ( 2.*out.row$rs.avg ) ) )
	out.row$rho.fit.a <-sqrt( 1. - exp( -2.*pmax( out.row$rs.avg - out.row$rs.fit, 0. ) / sum(expr.gains) ) )
	out.row$sum.rsq <- fit$sum.rsq
	out.row$sum.rsq0 <- fit.0$sum.rsq
	# put in interaction test
	out.row$rs.lin <- fit.lin$sum.res          # residual after linear model
	out.row$pv.asy <- va.asy$p.value              # p-value for asymmetries in model 

	 if (is.na(out.row$pv.fit)) {
	   stopifnot( fit.0$sum.rsq == 0. && fit$sum.rsq == 0. )
	   out.row$pv.fit <- 1.
	   out.row$xv.fit <- 0.
	 }
	
	 # append
	 outs[[i]] <- out.row

}


# concatenate outputs
source('toporbind.R')
out <- do.call( toporbind, outs )

# Normalized X2 value Data
range01 <- function(x, na.rm = T){(x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))}
normalized <- range01(out$xv.fit, na.rm = T)
out$xv.fit.01 <- normalized


# dump binary
saveRDS(out, 'model-stats-poi-nonlinear.rds')

# dump text
#write.gz.table <- function(data, filename, ...) {
#	file <- gzfile(filename, 'wb')
#	write.table(data, file, ...)
#	close(file)
#}
#write.gz.table(out, 'model-stats.tsv.gz', sep = '\t', row.names = T, col.names = NA)
