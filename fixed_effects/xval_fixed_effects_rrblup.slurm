#!/bin/bash
#SBATCH --job-name=xval_fixedeffects_rrblup # Name for your job
#SBATCH --account=uh
#SBATCH --partition=shared-long

#SBATCH --time=7-00:00:00 # Runtime in minutes. Default is 10 minutes. The Maximum runtime currently is 72 hours, 4320 minutes -requests over that time will not run

#SBATCH --nodes=1 # Number of nodes to spread cores across - default is 1 - if you are not using MPI this should likelybe 1
#SBATCH --tasks-per-node=1 # Number of tasks when using MPI. Default is 1
#SBATCH --cpus-per-task=1 # Number of cores requested, Default is 1 (total cores requested = tasks x cores)
#SBATCH --mem=100G # max amount of memory per node you require in MB (default) G (gigabyte)

#SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node

#SBATCH -e /home/quinnc3/Errors/%x-%A_%a.err # Standard err goes to this file
#SBATCH --mail-user quinnc3@hawaii.edu # this is the email you wish to be notified at
#SBATCH --mail-type ALL # this specifies what events you should get an email about ALL will alert you of jobbeginning,completion, failure etc
 
module load lang/R/4.4.1-gfbf-2023b

Rscript /home/quinnc3/fixed_effects/xval_fixed_effects_rrblup.R
