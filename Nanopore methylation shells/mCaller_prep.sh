#!/bin/bash
#SBATCH -J mCaller_prep
#SBATCH -p cpu
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --time=24:00:00
#SBATCH --mem=36G

### If you do not have access to these module replace with a single conda enviroment 
# Loading necessary modules and conda enviroment
module load miniconda/22.11.1-1
module load ont-guppy/6.5.7-cuda
module load bwa/0.7.17
module load gcc/11.2.0
module load samtools/1.14
conda activate nanopolish
###

# Default base name
base_name="Nanopore"

# Help function
function display_help() {
  echo "Usage: $0 [options] <input_FAST5_file> <input_fastq_file> <output_slow5_file>  
        <input_ref_file>"
  echo
  echo "   -h, --help                   Display this help message and exit."
  echo "   -n, --name <base_name>       Specify a base name for intermediate files. Default is 
          '$base_name'."
  echo "   <input_FAST5_file>           The input FAST5 file containing raw nanopore reads."
  echo "   <input_fastq_file>           The input FASTQ file, preferably processed through 
           Guppy for QA."
  echo "   <output_slow5_file>          The output file in SLOW5 format for downstream 
            processing."
  echo "   <input_ref_file>             The reference genome file in FASTA format for alignment 
            purposes."
  echo
  echo "This script processes nanopore sequencing data for mCaller analysis"
  exit 1
}

while [[ "$#" -gt 0 ]]; do
  case $1 in
    -h|--help) display_help ;;
    -n|--name) base_name="$2"; shift ;;
    --) shift; break ;;
    *) break ;;
  esac
  shift
done

# Check if the number of remaining arguments is less than expected
if [ "$#" -lt 4 ]; then
  echo "Error: Missing arguments."
  display_help
fi

# Assigning positional arguments after options
fast5_file="$1"
fastq_file="$2"
slow5_output="$3"
ref_genome="$4"

# Conversion from FAST5 to SLOW5 format
./slow5tools-v1.1.0/slow5tools f2s "$fast5_file" -o "$slow5_output"

# Indexing for Nanopolish
nanopolish index "$fast5_file" --slow5 "$slow5_output"

# Indexing the reference genome
bwa index "$ref_genome"

# Alignment and sorting
bwa mem -x ont2d -t 4 "$ref_genome" "$fastq_file" | samtools view -Sb - | samtools sort -T /tmp/"${base_name}".sorted -o "${base_name}".sorted.bam

samtools index "${base_name}".sorted.bam
# Event alignment
nanopolish eventalign -t 4 --scale-events -n -r "$fastq_file" -b "${base_name}".sorted.bam -g "$ref_genome" > "${base_name}".eventalign.tsv

# Ready for mCaller.sh 
echo "Processing completed. Intermediate files are named based on the base name: ${base_name}"



