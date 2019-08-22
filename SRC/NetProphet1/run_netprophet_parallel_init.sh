#!/bin/bash
#SBATCH --ntasks 4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH -D /scratch/mblab/levmorgan/NetProphet_2.0/SRC/NetProphet1
#SBATCH -o LOG/map_netprophet1_network_%A.out
#SBATCH -e LOG/map_netprophet1_network_%A.err

targetExpressionFile=${1}
regulatorExpressionFile=${2}
allowedMatrixFile=${3}
perturbationMatrixFile=${4}
differentialExpressionMatrixFile=${5}
microarrayFlag=${6}
nonGlobalShrinkageFlag=${7}
lassoAdjMtrFileName=${8}
combinedModelAdjMtrFileName=${9}
outputDirectory=${10}
combinedAdjLstFileName=${11}
regulatorGeneNamesFileName=${12}
targetGeneNamesFileName=${13}

echo "calling mpirun now, SLURM_NTASKS=${SLURM_NTASKS}"

mpirun -np ${SLURM_NTASKS} R --no-save -q --args ${targetExpressionFile} ${regulatorExpressionFile} ${allowedMatrixFile} ${perturbationMatrixFile} ${differentialExpressionMatrixFile} ${microarrayFlag} ${nonGlobalShrinkageFlag} ${lassoAdjMtrFileName} ${combinedModelAdjMtrFileName} ${outputDirectory} ${combinedAdjLstFileName} ${regulatorGeneNamesFileName} ${targetGeneNamesFileName} < run_netprophet_parallel_init.r
# mpirun -np ${NUM_JOBS} R --no-save -q --args ${targetExpressionFile} ${regulatorExpressionFile} ${allowedMatrixFile} ${perturbationMatrixFile} ${differentialExpressionMatrixFile} ${microarrayFlag} ${nonGlobalShrinkageFlag} ${lassoAdjMtrFileName} ${combinedModelAdjMtrFileName} ${outputDirectory} ${combinedAdjLstFileName} ${regulatorGeneNamesFileName} ${targetGeneNamesFileName} < run_netprophet_parallel_init.r
#srun --mpi=pmi2 Rscript /ufrc/data/training/SLURM/prime/rmpi_test.R
#srun --mpi=pmi2 R --no-save -q --args ${targetExpressionFile} ${regulatorExpressionFile} ${allowedMatrixFile} ${perturbationMatrixFile} ${differentialExpressionMatrixFile} ${microarrayFlag} ${nonGlobalShrinkageFlag} ${lassoAdjMtrFileName} ${combinedModelAdjMtrFileName} ${outputDirectory} ${combinedAdjLstFileName} ${regulatorGeneNamesFileName} ${targetGeneNamesFileName} < run_netprophet_parallel_init.r

