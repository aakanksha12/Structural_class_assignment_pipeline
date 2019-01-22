# Structural_class_assignment_pipeline

Structural class assignment pipeline assigns secondary structural classes and relative solvent accessibility classes to each site in a multiple sequence alignment using SCRATCH 1D structure prediction program.

## Getting Started
It has following requirements:
- SCRATCH1D prediction program (You need to provide the path for SCRATCH1D as a command line argument)
- A list of gene names for protein multiple sequence alignments as a text file
- Protein multiple sequence alignments in nexus format

## Usage

perl structure_assignment.pl gene_list.txt path_for_SCRATCH1D

## Output

It creates following output files
- A new nexus file with NEWMOD extension that has different Secondary structure and solvent accessibility classes site numbers written under a CHARSET block 
- SSPRO and ACCPRO output files 

### Important Note:

To extract sites for different structural classes from different genes (MSA), you can use PAUP to extract different CHARSET blocks from the output nexus file

### Sample files

Sample input and output files are provided in the example folder

