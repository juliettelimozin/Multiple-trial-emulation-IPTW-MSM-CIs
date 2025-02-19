#!/bin/bash
#SBATCH -J pseudo-true-300
#SBATCH -p icelake
#SBATCH --time=10:00:00
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=67
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk
#SBATCH --array=6,9,25

. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    

module load R/4.3.1-icelake

Rscript pseudo_true_300it_200k_parallel.R
