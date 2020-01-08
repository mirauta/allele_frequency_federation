# Allele_frequency_federation
This repository contains scripts that can be used to compute the allele frequencies using data generated in small cohorts.

# Download and running steps
The folowing instructions are for Windows users

## Install python with required packages
Download and install python - https://www.python.org/downloads/windows/

Download and install numpy - "pip install numpy" in the command line, ie CMD terminal

Download and install pandas - "pip install pandas" in the command line

## Download the  scripts

## Fill the samples_to_aggregate file
This should contain the folders and the names of the files on whiich partial MAF should be determined. If a file contains several samples, these need to be specified in the samples field (";" separated)
## Run
python genomes_to_maf.py --samples_file [samples_to_aggregate.xls] --destination [destination folder]
*replace [name] with the name of the file 

