#!/bin/bash
#SBATCH -A p30335                   # Allocation
#SBATCH -p long                     # Queue
#SBATCH -t 150:00:00                # Walltime/duration of the job
#SBATCH -N 1                        # Number of Nodes
#SBATCH --ntasks-per-node=6         # Number of Cores (Processors)
#SBATCH --mail-user=<jms@u.northwestern.edu  # Designate email address for job communications
#SBATCH --mail-type=END             # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --mail-type=FAIL            # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --output=sandbox/logs/sims.out           # Path for output must already exist
#SBATCH --error=sandbox/logs/sims.err            # Path for errors must already exist
#SBATCH --job-name="mrmi_sims"      # Name of job

echo "working directory = "$SLURM_SUBMIT_DIR

module load R/3.5.1
Rscript ./sandbox/simulations.R
