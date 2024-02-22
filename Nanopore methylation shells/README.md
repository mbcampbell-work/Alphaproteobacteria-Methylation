
This repository contains scripts used for the processing of raw nanopore data.

Note, that the latest verison of these programs are avaible at:
  for R.10+ data:       https://github.com/nanoporetech/dorado
  for R.9.4/9.5 data:   https://github.com/al-mcintyre/mCaller
    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dorado.sh
  Purpose: Converts raw nanopore sequencing data (Fast5) into methylation calls
  Notes: Only works with data post R.10

mCaller_prep.sh :
  Purpose: Pre-processing of raw r9.5 Nanopore sequencing data.
  Notes: Please run or use similar script before runing mCaller.sh

mCaller.sh:
  Purpose: Runs mCaller on provided aligned and sorted bam files
  Notes: Very long run time (atleast 24 hours, likely 48+ hours) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


