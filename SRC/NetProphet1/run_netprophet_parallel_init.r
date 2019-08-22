library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
targetExpressionFile <- toString(args[1])
regulatorExpressionFile <- toString(args[2])
allowedMatrixFile <- toString(args[3])
perturbationMatrixFile <- toString(args[4])
differentialExpressionMatrixFile <- toString(args[5])
microarrayFlag <- as.integer(args[6])
nonGlobalShrinkageFlag <- as.integer(args[7])
lassoAdjMtrFileName <- toString(args[8])
combinedAdjMtrFileName <- toString(args[9])
outputDirectory <- toString(args[10])
combinedAdjLstFileName <- toString(args[11])
regulatorGeneNamesFileName <- toString(args[12])
targetGeneNamesFileName <- toString(args[13])

source("run_netprophet_parallel.r")
library(Rmpi)
mpi.bcast.Robj2slave(args)

mpi.remote.exec(source("run_netprophet_parallel.r"))
# full.results <- mpi.remote.exec(run.netprophet.parallel(args))
full.results <- mpi.applyLB(1:8, function(job.num) {
															run.netprophet.parallel(list(job.num, args))
})
mpi.close.Rslaves()
save(full.results, file="full.results.raw.Rdata")
result.dfs <- list(bind_rows(lapply(full.results,
																		function(.res) {return(.res[[1]])})),
									 bind_rows(lapply(full.results,
																		function(.res) {return(.res[[2]])})))

result.dfs <- lapply(result.dfs, 
										 function(result) {
											 result.g <- group_by(result, regulators)
											 summarise_all(result.g, mean, na.rm=TRUE)
										 })

lasso_component <- subset(result.dfs[[1]], select=-c(regulators))
write.table(lasso_component,file.path(outputDirectory,lassoAdjMtrFileName),row.names=FALSE,col.names=FALSE,quote=FALSE)
combinedAdjMtr <- subset(result.dfs[[2]], select=-c(regulators))
write.table(combinedAdjMtr,combinedAdjMtrFileName,row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

cat("Successfully generated solution at ", outputDirectory, "\n")
mpi.quit()
#full.results.df <- bind_rows(full.results)



## Perform model averaging to get final NetProphet Predictions

# if(length(args) == 13 & file.exists(regulatorGeneNamesFileName) & file.exists(targetGeneNamesFileName)){
#	 source("make_adjacency_list.r")
# }

# Bootstrapping approach:
# We have to run model selection 20x - write a separate R file with relevant logic
# for each beta in each model, evaluate s_i,j = \sigma^2 without that beta (set it to 0?)/\sigma^2 with it
# average score for each beta, this is our final ranked list


# Use mpi.applyLB over a list/vector of genes
# function takes a gene name and whatever data we're running lasso on
# runs boot to calculate variance explained with the gene and without the gene 20x
