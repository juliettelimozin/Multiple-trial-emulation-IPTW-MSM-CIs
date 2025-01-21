#!/bin/bash
#SBATCH -J true-boot
#SBATCH -p icelake
#SBATCH --time=10:00:00
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk

. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    

module load R/4.3.1-icelake

Rscript true_value_MSM_bootstrap_200it.R
