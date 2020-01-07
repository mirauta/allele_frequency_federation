# allele_frequency_federation
This repository contains scripts that will enable the federation of allele frequencies using data generated in small cohorts.

## Download the  scripts
## Install required packages
python - https://www.python.org/downloads/windows/

numpy - "pip install numpy" in the command line, ie CMD terminal

pandas - "pip install pandas" in the command line

## Fill the samples_to_aggregate file
This should contain the folders and the names of the files on whiich partial MAF should be determined. If a file contains several samples, these need to be specified in the samples field (";" separated)
## Run
python genomes_to_maf.py --samples_file samples_to_aggregate.xls --destination [destination folder]

