#!/bin/bash -l

#SBATCH -A snic2017-7-16                        # Uppmax project
#SBATCH -p core                                 # choice between cores or nodes
#SBATCH -n 8                                   # number of cores
#SBATCH -t 8:00:00                              # computing time
#SBATCH -J mkfastq_job                          # job name
#SBATCH --mail-user panosbino@gmail.com
#SBATCH --mail-type=ALL   


screen

#load required modules
module load bioinfo-tools
module load cellranger

# generate FASTQs from BCL files found in Illumina Run folder

cellranger mkfastq --id=10x95 --run=/proj/uppstore2018019/private/raw/210414_NB502120_0265_AH3VHMBGXJ --csv=/proj/uppstore2018019/private/raw/210414_NB502120_0265_AH3VHMBGXJ/input_samplesheet2.csv


