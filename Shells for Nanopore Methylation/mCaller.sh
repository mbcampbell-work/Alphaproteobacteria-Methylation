#!/bin/bash
#SBATCH -J mCaller
#SBATCH -p gpu
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH -G 2 # Number of GPUS 
#SBATCH -c 4  # Number of Cores per Task
#SBATCH --time=18:00:00
#SBATCH --mem=80G

#########################################
# RUN mCaller_prep.sh before running	#
#########################################

module load miniconda/22.11.1-1
module load ont-guppy/6.5.7-cuda
module load bwa/0.7.17
module load gcc/11.2.0
module load samtools/1.14

conda activate mcaller2

# Replace <file> with the name chosen during the prep stage
# Use same reference that was indexed in prep stage
python3 mCaller/mCaller.py \
  -m A \
  -b A \
  -r <reference>.fasta \
  -e <file>.eventalign.tsv \
  -d mCaller/r95_twobase_model_NN_6_m6A.pkl \
  -f <file>.fastq 
  -t 2
 
 python ./mCaller/make_bed.py -f <file>.eventalign.diffs.6 -d 40 -t 0 --ref <reference>.fasta --vo