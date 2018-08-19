# Structural_class_assignment_pipeline

Structural class assignment pipeline assigns secondary sructural classes and relative solvent accessibility classes to each site in a multiple sequence alignment using SCRATCH 1D structure prediction program.

## Getting Started
It has following requirements:
- SCRATCH1D prdiction program (You have to provide the path for SCRATCH1D as a command line argument)
- A list of gene names for protein multiple sequence alignments as a text file
- Protein sequence alignments in nexus format

## Usage

perl structure_assign.pl gene_list.txt current_working_directory path_for_SCRATCH1D

## Output

It creates following output files
- A new nexus file with NEWMOD extension that has different Secondary structure and solvent accessibility classes site numbers written under a CHARSET block 
- SSPRO and ACCPRO output files 

### Important Note:

In order to extract different structural classes sites from different genes, you need PAUP to extract differnt CHARSET blocks from the output nexus file

### Sample file

Sample input and output files are provided in the example folder

