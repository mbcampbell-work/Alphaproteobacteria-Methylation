# IGV_Brucella Data

This folder contains data and instructions for visualizing Brucella genome data using the Integrative Genomics Viewer (IGV) desktop application. Please note that the Quick Start process will not work in the web browser version of IGV.

## Quick Start for comparing strains

1. Load IGV desktop.
2. Click on `File` in the top toolbar.
3. Select `Open Session`.
4. Navigate to the directory where you have saved the IGV_Brucella Files.
5. Click on the file named `igv_comparison.xml`.

## Definitions

- [Strain] VS WT wig:
   - Shows the probability of methylation at each GANTC when comparing the probability of the wild type vs the strain. Positive values (blue) signifies that the strain in question has a higher probability of methylation compared to wild type, while negative values (red) signifies that wild type has a higher probability of methylation compared to the other strain.

- MUCRbinding.wig:
   - Peaks show regions that are considered binding sites for MucR.

- MUCRbindingsites.bed:
   - Bands show regions that are considered binding sites for MucR.

## Quick Start for inspecting individual strains.

1. Load IGV desktop.
2. Click on `File` in the top toolbar.
3. Select `Open Session`.
4. Navigate to the directory where you have saved the IGV_Brucella Files.
5. Click on the file named `igv_siite_by_site.xml`.

Definitions

- [Strain].wig:
   - Shows the probability of methylation at each GANTC site, positive values(blue) signify GANTC sites on the plus strand and negative values (red) signify GANTC on the negative strand. 

## Loading File Outside of `igv.xml`

1. Load IGV (any version).
2. Click on `Genome` in the top toolbar.
3. Select `Load Genome from Files`.
4. Navigate to the directory where you have saved the IGV_Brucella Files.
5. Select all files within the Genome subfolder at the same time and select confirm (`ref.fasta` and `ref.fasta.fai`)
6. Additional information can be loaded by files from each of the subfolders below:

## Subfolders

The IGV_Brucella data is organized into 4 subfolders. In order to load anything outside of the `igv.xml` file, you must load all of the files from the `Genome Fasta` folder (`ref.fasta` and `ref.fasta.fai`).


1. Genome Fasta:
Contains ref.fasta and res.fasta.fai files required for load genome option. Note that both must be uploaded together. 

2. Gene Features
Contains three files: Chromosome 1.gff3, Chromosome 2.gff3, mucRbindingistes.wig . Loading either of the gff3 files will allow you to search up genes with either BAB_1 format, BAB_RS format, or gene name. Mucrbindingsites.wig will add bands to regions defined as MucR binding sites.

3. Site by Site:
Contains four wig files, that when loaded will display the site by site probability of methylation for each strain. 

4. Comparing Strains:
Contains three wig files, each of which shows the probability of methylation at each GANTC when comparing the probability of the wild type vs the strain. Positive values (blue) signifies that the strain in question has a higher probability of methylation compared to wild type, while negative values (red) signifies that wild type has a higher probability of methylation compared to the other strain.

