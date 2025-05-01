#!/bin/bash
#SBATCH -J simu_high_small
#SBATCH -p icelake
#SBATCH --time=3:0:0
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=67
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk
#SBATCH --array=1-9


. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    


module load R/4.3.1-icelake


Rscript simulation_high_jackknife.R
