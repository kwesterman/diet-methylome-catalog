#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --output=slurm-numbers.%N.%j.$1.out 
#SBATCH --error=slurm-numbers.%N.%j.err
#SBATCH --mail-user=$USER

if [[ $# -ne 2 ]]	# Check for exactly two arguments
	then
		echo "Must supply two arguments: exposure_var and cov_suffix"
		exit 1
fi

module load R/3.3.2 	# Load R module
Rscript --no-save gem_emodel.R $1 $2	# Call R script to format data and run GEM analysis
