#!/bin/bash
# RCTD analysis with new R environment

# Partition for the job:
#SBATCH --partition=long

# Multithreaded (SMP) job: must run on one node
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="rctd_liver_newenv"

# The project ID which this job should run under:
#SBATCH --account="punim1236"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=100G

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=31-0:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from this directory:
cd /data/projects/punim1236/liver_stomicsv1.3/2026714RCTDnew

# Load R 4.5.0 with new library
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0
export R_LIBS_USER=/data/projects/punim1236/liver_stomicsv1.3/RCTD_CLAHE_old/R_library_4.5

# The job command:
Rscript /data/projects/punim1236/liver_stomicsv1.3/2026714RCTDnew/RCTD_CLAHE_bin20_liverSTOmics_newenv.R

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s
