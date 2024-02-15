#!/bin/bash
#SBATCH -J Dorado
#SBATCH -p gpu # Note requires gpu nodes
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH -G 4 # Number of GPUS 
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --time=24:00:00
#SBATCH --mem=60G

module load ont-guppy/6.5.7-cuda
module load miniconda/22.11.1-1
module load openmpi/4.1.4+cuda11.6.2-ucx
module load dorado/0.3.1+cuda11.6.2
module load samtools/1.9

# DO NOT CHANGE, unless using new nanopore data (dependent on pore) 
BASECALLING_MODEL="dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
MODBASES_MODEL="dna_r10.4.1_e8.2_400bps_sup@v4.2.0_6mA@v2"

# Change to where you have your input files
RUN_NAME="FINAL_RUN"
# Requires references fasta file of genome in RUN_NAME file
REFERENCE="ref.fasta"

mkdir "$RUN_NAME/input"

# Need to either pip install or create pod5 enviroment 
# https://pod5-file-format.readthedocs.io/en/latest/docs/install.html
conda activate pod5 

pod5 convert fast5 ../Nanopore/QUO1007176/NanoporeDNAReads/FAST5/2308/*.fast5 --output "$RUN_NAME/input/WT_Brucella/output.pod5"

pod5 convert fast5 ../Nanopore/QUO1007176/NanoporeDNAReads/FAST5/2308pRW414/*.fast5 --output "$RUN_NAME/input/2308pRW414_Brucella/output.pod5"

pod5 convert fast5 ../Nanopore/QUO1007176/NanoporeDNAReads/FAST5/GR106/*.fast5 --output "$RUN_NAME/input/GR106_Brucella/output.pod5"

pod5 convert fast5 ../Nanopore/QUO1007176/NanoporeDNAReads/FAST5/CC092/*.fast5 --output "$RUN_NAME/input/CC092_Brucella/output.pod5"

mkdir "$RUN_NAME/intermediate"
mkdir "$RUN_NAME/output"

POD_FILES=("$RUN_NAME/input/WT_Brucella/" "$RUN_NAME/input/GR106_Brucella/" "$RUN_NAME/input/CC092_Brucella/" "$RUN_NAME/input/2308pRW414_Brucella/"
)
BAM_FILES=("$RUN_NAME/intermediate/WT_Brucella.bam" "$RUN_NAME/intermediate/GR106_Brucella.bam" "$RUN_NAME/intermediate/CC092_Brucella.bam" "$RUN_NAME/intermediate/2308pRW414_Brucella.bam")

BAM_NAMES=("WT_Brucella.bam" "GR106_Brucella.bam" "CC092_Brucella.bam" "2308pRW414_Brucella.bam")

SixmA_BED_NAMES=("6mA_WT_Brucella.bed" "6mA_GR106_Brucella.bed" "6mA_CC092_Brucella.bed" "6mA_2308pRW414_Brucella.bed")
    
for ((i=0; i<${#POD_FILES[@]}; i++)); do
    dorado-0.3.4-linux-x64/bin/dorado basecaller  \
    $BASECALLING_MODEL \
    "${POD_FILES[$i]}" \
    --modified-bases-models $MODBASES_MODEL \
    --device cuda:all \
    --verbose \
    > "${BAM_FILES[$i]}"
done

for ((i=0; i<${#BAM_FILES[@]}; i++)); do
    dorado aligner $REFERENCE ${BAM_FILES[$i]} >  "$RUN_NAME/intermediate/aligned_${BAM_NAMES[$i]}"
    echo "Done aligning ${BAM_NAMES[$i]}"
done

for ((i=0; i<${#BAM_FILES[@]}; i++)); do
    samtools sort "$RUN_NAME/intermediate/aligned_${BAM_NAMES[$i]}" -o "$RUN_NAME/intermediate/sorted_${BAM_NAMES[$i]}"
    echo "Done sorting ${BAM_NAMES[$i]}"
done

for ((i=0; i<${#BAM_FILES[@]}; i++)); do
    samtools index "$RUN_NAME/intermediate/sorted_${BAM_NAMES[$i]}"
    echo "Done indexing ${BAM_NAMES[$i]}"
done

# Need to download and create modbam2bed enviroment
# Note that modkit also works 
conda activate modbam2bed

for ((i=0; i<${#BAM_FILES[@]}; i++)); do

    modbam2bed --mod_base=6mA $REFERENCE "$RUN_NAME/intermediate/sorted_${BAM_NAMES[$i]}" > "$RUN_NAME/output/${SixmA_BED_NAMES[$i]}"
    
    echo "Done creating ${SixmA_BED_NAMES[$i]}"
    
done


