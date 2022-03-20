# Sim Genomes

Given a FASTA file, the script simulates the specified number of SNPs and indels within the genome. The simulated variants are also output in a VCF file. 
    
&nbsp;

&nbsp;


## Usage

Takes the no. of SNPs to be induced, the no. of indels to be induced, the reference FASTA file and the name of the output file as input arguments. Outputs the genome sequence of the simulated sequence in a FASTA file and the variants simulated in a VCF file.

The script should be run with the following commands:
```python
python vcfeval_filter.py snps indels seq_file out_file
```
