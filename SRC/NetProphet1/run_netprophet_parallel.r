run.netprophet.parallel <- function(args) {
	my.rank <- mpi.comm.rank()
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


	cat("Loading data...\n")
	gene.names <- make.names(read.table(targetGeneNamesFileName,
													 stringsAsFactors=F)[,1])
	reg.names <- make.names(read.table(regulatorGeneNamesFileName,
													stringsAsFactors=F)[,1])

	# Target data. Rows are genes, columns are samples
	tdata <- as.matrix(read.table(targetExpressionFile, row.names=gene.names))
	# Regulator data. Rows are TFs, columns are samples
	rdata <- as.matrix(read.table(regulatorExpressionFile, row.names=reg.names))
	allowed <- as.matrix(read.table(allowedMatrixFile, row.names=reg.names,
																	col.names=gene.names))
	pert <- as.matrix(read.table(perturbationMatrixFile, row.names=gene.names))
	de_component <- as.matrix(read.table(differentialExpressionMatrixFile,
																			 row.names=reg.names,
																			 col.names=gene.names))

	source("run_netprophet.r")

	res <- run_netprophet(gene.names,
												reg.names,
												tdata,
												rdata,
												allowed,
												pert,
												de_component)
	lasso_component <- res[[1]]
	write.table(lasso_component,file.path(outputDirectory,paste0("lasso.",rank,".Rdata")),row.names=FALSE,col.names=FALSE,quote=FALSE)
	combinedAdjMtr <- res[[2]]
	write.table(combinedAdjMtr,file.path(outputDirectory,paste0("combined.",rank,".Rdata")),row.names=FALSE,col.names=FALSE,quote=FALSE)

	cat("Successfully generated solution at ", outputDirectory, "\n")
	return(res)
}
