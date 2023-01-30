#!/bin/bash
#SBATCH -J size-red
#SBATCH -p icelake
#SBATCH --time=3:0:0
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=67
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk
#SBATCH --array=2,5,8,10,25,50

. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    

module load R/4.2.0-icelake

Rscript est_size_red.R

