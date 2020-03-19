#!/bin/bash
#
# Re-run the code locally, to re-create the data and figure.
#
# Usage:
#
#   ./scripts/rerun.sh
#
#SBATCH --partition=gelifes
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --job-name=pirex24
#SBATCH --output=example_24.log
#
rm -rf example_24
rm *.png
time Rscript example_24.R
zip -r pirouette_example_24.zip example_24 example_24.R scripts *.png

