#!/bin/bash
#SBATCH -J jackknife_high_fixed_large
#SBATCH -p icelake
#SBATCH --time=17:0:0
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=67
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk
#SBATCH --array=19


. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    


module load R/4.3.1-icelake


Rscript simulation_high_jackknife.R
