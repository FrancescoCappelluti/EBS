#!/bin/bash
#SBATCH --account=
#SBATCH --nodes=4
#SBATCH --mem=
#SBATCH --ntasks-per-node=
#SBATCH --cpus-per-task=1
#SBATCH --job-name=orca_EBS
#SBATCH --time=24:00:00

module load mkl
INPUT=Fe4S4
OUTPUT=$INPUT
OUTPUT_ERR=out.err

ROOT=$SLURM_SUBMIT_DIR

cd $ROOT
scontrol show hostnames $SLURM_JOB_NODELIST > machinefile

echo "Starting run at: `date`"
echo "Job ID:          $SLURM_JOBID"
echo "Job name:        $SLURM_JOB_NAME"
echo "Output file:     $ROOT/$OUTPUT"
echo "PBS Node file(machine file):     $PBS_NODEFILE"

cd ${ROOT}/geo_opt/

orca ${INPUT}.inp > ${OUTPUT}.out 2> ${OUTPUT_ERR} 
