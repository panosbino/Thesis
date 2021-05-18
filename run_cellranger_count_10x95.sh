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


cellranger count --id=10x95 --transcriptome=/proj/uppstore2018019/private/panos/Mouse_genome_pMR641 --fastqs=/proj/uppstore2018019/private/raw/210414_NB502120_0265_AH3VHMBGXJ/10x95/outs/fastq_path/H3VHMBGXJ/10x95 --sample=10x95 