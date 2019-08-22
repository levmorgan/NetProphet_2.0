run.netprophet.parallel <- function(job) {
	my.rank <- job[[1]]
	args <- job[[2]]
	set.seed(747+my.rank)
	targetExpressionFile <- toString(args[1])
	regulatorExpressionFile <- toString(args[2])
	allowedMatrixFile <- toString(args[3])
	perturbationMatrixFile <- toString(args[4])
	differentialExpressionMatrixFile <- toString(args[5])
	microarrayFlag <- as.integer(args[6])
	nonGlobalShrinkageFlag <- 1# as.integer(args[7])
	lassoAdjMtrFileName <- toString(args[8])
	combinedAdjMtrFileName <- toString(args[9])
	outputDirectory <- toString(args[10])
	combinedAdjLstFileName <- toString(args[11])
	regulatorGeneNamesFileName <- toString(args[12])
	targetGeneNamesFileName <- toString(args[13])

	source("run_netprophet.r")

	cat("Loading data...\n")
	# Create fake regulator and gene names for BBSR
	# We're not using real ones because the data may be bootstrapped, 
	# resulting in collisions
	real.reg.names <- make.names(read.table(regulatorGeneNamesFileName,
													stringsAsFactors=F)[,1])
	real.gene.names <- make.names(read.table(targetGeneNamesFileName,
													 stringsAsFactors=F)[,1])
	reg.names <- paste0("reg.", 1:length(real.reg.names))
	gene.names <- paste0("gene.", 1:length(real.gene.names))


	# Target data. Rows are genes, columns are samples
	tdata <- as.matrix(read.table(targetExpressionFile, row.names=gene.names))
	# Regulator data. Rows are TFs, columns are samples
	rdata <- as.matrix(read.table(regulatorExpressionFile, row.names=reg.names))

	# Resample regulator data for bootstrapping
	regulator.idx <- sample(1:dim(rdata)[1], replace=T)
	rdata <- rdata[regulator.idx,]

	allowed <- as.matrix(read.table(allowedMatrixFile, row.names=reg.names,
																	col.names=gene.names))
	pert <- as.matrix(read.table(perturbationMatrixFile, row.names=gene.names))
	de_component <- as.matrix(read.table(differentialExpressionMatrixFile,
																			 row.names=reg.names,
																			 col.names=gene.names))


	res <- run_netprophet(reg.names,
												gene.names,
												tdata,
												rdata,
												allowed,
												pert,
												de_component,
												microarrayFlag,
												bootstrap.run=T
												)
	lasso_component <- res[[1]]
	write.table(lasso_component,file.path(outputDirectory,paste0("lasso.",my.rank,".csv")),row.names=FALSE,col.names=FALSE,quote=FALSE)
	combinedAdjMtr <- res[[2]]
	write.table(combinedAdjMtr,file.path(outputDirectory,paste0("combined.",my.rank,".csv")),row.names=FALSE,col.names=FALSE,quote=FALSE)
	res[[1]] <- as.data.frame(res[[1]])
	res[[2]] <- as.data.frame(res[[2]])
	res[[1]]$regulators <- reg.names
	res[[2]]$regulators <- reg.names

	cat("Successfully generated partial solution at ", outputDirectory, "\n")
	return(res)
}
