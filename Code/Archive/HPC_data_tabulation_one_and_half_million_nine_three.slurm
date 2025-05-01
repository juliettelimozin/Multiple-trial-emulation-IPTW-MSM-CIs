#!/bin/bash
#SBATCH -J true-boot-test
#SBATCH -p icelake
#SBATCH --time=1:00:00
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk

. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    

module load R/4.3.1-icelake

Rscript data_tabulation_one_and_half_million_nine_three.R
