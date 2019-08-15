library(future.apply)
plan(multiprocess, workers=availableCores()-1)

make.matrix.mask <- function(v1, v2) {
	mask <- as.logical(outer(v1, v2))
	mask <- matrix(mask, nrow=length(v1), ncol=length(v2))
	mask
}

run_netprophet <- function(gene.names,
													 reg.names,
													 tdata,
													 rdata,
													 allowed,
													 pert,
													 de_component) {
	targets <- seq(dim(tdata)[1])

	if(microarrayFlag == 0) {
		##RNA-Seq Data
		tdata <- log(tdata+1)/log(2)
		rdata <- log(rdata+1)/log(2)
	}

	## Skip regression on some genes
	cat("Gene samples skipped count:\n")
	skip_gen <- rep(0, dim(tdata)[1])
	for (i in 1:dim(tdata)[1]) {
		if (sum(tdata[i,] != 0) < dim(tdata)[2]/10+1) {
			skip_gen[i] = 1
		}
	}
	# skip_gen <- future_apply(tdata, 1, function(row) {
	# 	if (sum(row != 0) < dim(tdata)[2]/10+1) {
	# 		return(1)
	# 	}
	# 	return(0)
	# })
	cat(length(which(skip_gen == 1)), "\nRegulator samples skipped:\n")
	skip_reg <- rep(0, dim(rdata)[1])
	for (i in 1:dim(rdata)[1]) {
		if (sum(rdata[i,] != 0) < 1) {
			skip_reg[i] = 1
			cat(i, ",")
		}
	}
	cat("\n")

	cat("Centering and scaling data\n")
	## Center data
	tdata <- tdata - apply(tdata,1,mean)
	rdata <- rdata - apply(rdata,1,mean)

	## Scale data
	t.sd <- apply(tdata,1,sd)
	t.sdfloor <- mean(t.sd) + sd(t.sd)
	t.norm <- apply(rbind(rep(t.sdfloor,times=length(t.sd)),t.sd),2,max) / sqrt(dim(tdata)[2])
	tdata <- tdata / ( t.norm * sqrt(dim(tdata)[2]-1) )
	# tdata <- tdata / t.sd
	#
	r.sd <- apply(rdata,1,sd)
	r.sdfloor <- mean(r.sd) + sd(r.sd)
	r.norm <- apply(rbind(rep(r.sdfloor,times=length(r.sd)),r.sd),2,max) / sqrt(dim(rdata)[2])
	rdata <- rdata / ( r.norm * sqrt(dim(rdata)[2]-1) )
	# rdata <- rdata / r.sd

	## Compute unweighted solution
	cat("Computing unweighted solution\n")
	prior <- matrix(1,ncol=dim(tdata)[1], nrow=dim(rdata)[1],
									dimnames=list(reg.names, gene.names))

	## TODO: Both seed and # of cv folds (in global.lars.regulators.r) are parameters that should be exposed to the user
	seed <- 747
	set.seed(seed)

	cat("Computing mutual information between regulators and targets\n")
	x <- rdata * prior[,i]

	# Remove disallowed genes and regulators
	filt.tdata <- t(tdata[,skip_gen == 0])
	filt.rdata <- t(rdata[,skip_reg == 0])

	source("mixed_clr.r")
	mutual.information <- t(mi(x=filt.tdata, y=filt.rdata, cpu.n=availableCores()))
	cat("Computing CLR\n")
	clr.results <- mixedCLR(mi.stat=data.frame(), mi.dyn=mutual.information)
	clr.results <- t(clr.results)
	rownames(clr.results) <- gene.names
	colnames(clr.results) <- reg.names
	clr.rankings <- apply(clr.results, 2, function(col) {order(col, decreasing=T)})

	source("bayesianRegression.R")
	cat("Computing Bayesian Best Subset Regression\n")
	models <- BBSR(X=rdata, Y=tdata, clr.mat=clr.results,
								 nS=10, no.pr.val=0, weights.mat=t(allowed),# Come back to this
								 prior.mat=t(prior), cores=availableCores()-1)
	# bbsr.results has regulators as rows and genes as columns
	bbsr.results <- matrix(0,ncol=dim(tdata)[1], nrow=dim(rdata)[1],
									dimnames=list(reg.names, gene.names))
	for (i in 1:length(models)) {
		gene.model <- models[[i]]
		bbsr.results[,i][gene.model$pp] <- gene.model$betas.resc
	}
		
	lasso_component <- bbsr.results

	# write.table(lasso_component,file.path(outputDirectory,lassoAdjMtrFileName),row.names=FALSE,col.names=FALSE,quote=FALSE)

	## Perform model averaging to get final NetProphet Predictions
	source("combine_models.r")
	combinedAdjMtr <- combine_models(lasso_component)
	return(list(lasso_component, combinedAdjMtr)) 
}
