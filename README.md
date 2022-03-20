# SeqIdCalc

Reclassifies the variants output by vcfeval that been incorrectly evaluated. It also checks for incorrect evaluation due to left-alignment of variants thus requires the reference FASTA file. New VCF files are written containing the corrected true positive, false positive and false negative calls. 
    
&nbsp;

&nbsp;


## Usage

Takes the reference FASTA file, truthset VCF file, calls VCF file and vcfeval directory path as input. Outputs the corrected TP, FP and FN VCF files within the vcfeval directory.

The script should be run with the following commands:
```python
python vcfeval_filter.py reference truthset called_vcf vcfeval_dir
```
