#!/bin/bash

#SBATCH --job-name=05_Build_SDMs
#SBATCH --output=/work/%u/%x-%A-%a.log
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-type=BEGIN,END,FAIL

module load foss/2020b R/4.0.4-2

taxon_name="$1"
species_csv="data_raw/$taxon_name.csv"

output_dir="/work/$USER/$SLURM_JOB_NAME-$SLURM_JOB_ID"
mkdir -p "$output_dir"


Rscript src/05_Build_SDMs.R "$taxon_name" "$species_csv" "$output_dir"

# cd ~/SoilBiodiversity
# sbatch -a 1-$(xsv count data_raw/Crassiclitellata.csv) src/05_Build_SDMs.sh Crassiclitellata 



