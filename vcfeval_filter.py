'''
3/5/21
Python script to reclassify the variants output by vcfeval that been incorrectly evaluated. Takes the reference FASTA file, truthset VCF file, calls VCF file and vcfeval directory path as input. Also, checks for incorrect evaluation due to left-alignment of variants. Outputs the corrected TP, FP and FN VCF files within the vcfeval directory.  
'''

import os
import sys
#import vcf
import gzip
from itertools import chain
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('reference', type=str, 
                   help='Reference FASTA file')
parser.add_argument('truthset', type=str, 
                   help='Truthset VCF file')
parser.add_argument('called_vcf', type=str, 
                   help='Called VCF file')
parser.add_argument('vcfeval_dir', type=str, 
                   help='vcfeval dir')

args = parser.parse_args()


fp_file = 'fp.vcf.gz'
new_tp = 'tp_corrected.vcf.gz'
old_tp = 'tp.vcf.gz'
new_fp = 'fp_corrected.vcf.gz'
old_fn = 'fn.vcf.gz'
new_fn = 'fn_corrected.vcf.gz'


#read in truthset
truthset_path = args.truthset
#truthset_broken = []
truthset_broken = [l.decode("utf-8") for l in gzip.open(truthset_path, "rb").readlines() if not l.startswith(b'#')]

#create an array to contain the full line in vcf
truthset_full = [None] * len(truthset_broken)

#read in lines and remove additional info
for index, line in enumerate(truthset_broken):
    line_list = list(line.split('\t'))
    #write for line for later editing into truthset_full
    truthset_full[index] = line_list
    new_list = []
    new_list.append(line_list[1])
    new_list.append(line_list[3])
    new_list.append(line_list[4])
    truthset_broken[index] = new_list


#fn file
fn_path = os.path.join(args.vcfeval_dir, old_fn)
new_fn_path = os.path.join(args.vcfeval_dir, new_fn)
fn_cor = gzip.open(new_fn_path, "w")
fn_cor.writelines([l for l in gzip.open(fn_path, "r").readlines() if l.startswith(b"#")])

#fp file
fp_path = os.path.join(args.vcfeval_dir, fp_file)
#specify new file path
new_fp_path = os.path.join(args.vcfeval_dir, new_fp)
fp_cor = gzip.open(new_fp_path, "w")
fp_cor.writelines([l for l in gzip.open(fp_path, "r").readlines() if l.startswith(b"#")])

#create a new tp file
new_tp_path = os.path.join(args.vcfeval_dir, new_tp)
#write old tp file to corrected tp file
old_tp_path = os.path.join(args.vcfeval_dir, old_tp)
#write the tp file into the new file
tp_cor = gzip.open(new_tp_path, "w")
tp_cor.writelines([l for l in gzip.open(old_tp_path, "r").readlines()])
tp_read = [[l.decode("utf-8").split('\t')[1],l.decode("utf-8").split('\t')[3], l.decode("utf-8").split('\t')[4]] for l in gzip.open(old_tp_path, "r").readlines() if not l.startswith(b"#")]


#read in the reference FASTA file
with open(args.reference, "r") as r:
    #reads in the sequence
    ref_seq = list("".join(l.strip() for l in r.readlines() if not str(l).startswith('>')))



fp_file_list = []

#opens fp file and truthset for reading 
with gzip.open(fp_path,"r") as f:
    for line in f:
        #writes out header to new fp file
        if not line.startswith(b'#'):
            fp_line = list(line.decode("utf-8").split('\t'))

            #checks if variant position isn't in the truthset 
            if not {fp_line[1]}.issubset(chain.from_iterable(truthset_broken)):
                #checks insertions
                if len(fp_line[4]) >= len(fp_line[3]): 
                    #finds possible variants that may be represented differently through left-alignment
                    variant_start = int(fp_line[1]) - (2*len(fp_line[4]))
                    variant_end = int(fp_line[1]) + (2*len(fp_line[4]))

                    #checks for possible positions of variant in truthset
                    for var in range(variant_start, variant_end+1):
                        left_aligned = [t for t in truthset_broken if (int(t[0]) == int(var)) and (int(len(t[2])) > int(len(t[1])))]
                        if left_aligned:
                            break

                    #if there is no match in truth set, writes out variant as false positive 
                    if not left_aligned:
                        fp_cor.write(line)
                        fp_file_list.append(line)

                    #checks if the variants in truthset and FP vcf are the same
                    else:  
                        seq_end_pos = int(max(fp_line[1], left_aligned[0][0])) + len(fp_line[4])

                        #get the original sequence
                        sequence = ref_seq[int(min(fp_line[1], left_aligned[0][0]))-1:seq_end_pos]

                        #recreate insertion from the FP VCF file and truthset to compare if they are the same call
                        #compares positions to see if the FP call is more left-aligned
                        if int(fp_line[1]) < int(left_aligned[0][0]):
                            fp_sequence = list(str(fp_line[4])) + sequence[len(fp_line[3]):]

                            pos_diff = int(left_aligned[0][0]) - int(fp_line[1])
                            truth_sequence = sequence[:pos_diff] + list(str(left_aligned[0][2])) + sequence[pos_diff+len(left_aligned[0][1]):]

                            #if the variants are the same, writes out to  
                            if fp_sequence == truth_sequence:
                                tp_cor.write(line)
                                tp_read.append([fp_line[1], fp_line[3], fp_line[4]])
                                tp_read.append(left_aligned[0])
                            else:
                                fp_cor.write(line)
                                fp_file_list.append(line) 

                        #compares positions to see if the truthset is more left-aligned
                        if int(fp_line[1]) > int(left_aligned[0][0]):
                            truth_sequence = list(str(left_aligned[0][0])) + sequence[len(left_aligned[0][2]):]

                            pos_diff = int(fp_line[1]) - int(left_aligned[0][0])
                            fp_sequence = sequence[:pos_diff] + list(str(fp_line[4])) + sequence[pos_diff+len(fp_line[3]):]


                            #if the variants are the same, writes out to  
                            if fp_sequence == truth_sequence:
                                    tp_cor.write(line)
                                    tp_read.append([fp_line[1], fp_line[3], fp_line[4]])
                                    tp_read.append(left_aligned[0])
                            else:
                                fp_cor.write(line)
                                fp_file_list.append(line)  


                #checks deletions
                if len(fp_line[3]) > len(fp_line[4]): 
                    #finds possible variants that may be represented differently through left-alignment
                    variant_start = int(fp_line[1]) - (2*len(fp_line[3]))
                    variant_end = int(fp_line[1]) + (2*len(fp_line[3]))

                    #checks for possible positions of variant in truthset
                    for var in range(variant_start, variant_end+1):
                        left_aligned = [t for t in truthset_broken if (int(t[0]) == int(var)) and (int(len(t[1])) > int(len(t[2])))]
                        if left_aligned:
                            break

                    #if there is no match in truth set, writes out variant as false positive 
                    if not left_aligned:
                        fp_cor.write(line)
                        fp_file_list.append(line)

                    #checks if the variants in truthset and FP vcf are the same
                    else: 
                        seq_end_pos = int(max(fp_line[1], left_aligned[0][0])) + max(len(fp_line[3]), len(left_aligned[0][1]))

                       #get the original sequence
                        sequence = ref_seq[int(min(fp_line[1], left_aligned[0][0]))-1:seq_end_pos]

                        #recreate insertion from the FP VCF file and truthset to compare if they are the same call
                        if int(fp_line[1]) < int(left_aligned[0][0]):
                            fp_sequence = sequence[:len(fp_line[4])] + sequence[len(fp_line[3]):]

                            pos_diff = int(left_aligned[0][0]) - int(fp_line[1])
                            truth_sequence = sequence[:(pos_diff+len(left_aligned[0][2]))] + sequence[(pos_diff+len(left_aligned[0][1])):]

                            #if the variants are the same, writes out to TP vcf file  
                            if fp_sequence == truth_sequence:
                                tp_cor.write(line)
                                tp_read.append([fp_line[1], fp_line[3], fp_line[4]])
                                tp_read.append(left_aligned[0])

                            else:
                                fp_cor.write(line)
                                fp_file_list.append(line)

                        #checks if FP call position is greater than truthset position
                        if int(fp_line[1]) > int(left_aligned[0][0]):
                            truth_sequence = sequence[:len(left_aligned[0][2])] + sequence[len(left_aligned[0][1]):]

                            pos_diff = int(left_aligned[0][0]) - int(fp_line[1])
                            fp_sequence = sequence[:(pos_diff+len(fp_line[4]))] + sequence[(pos_diff+len(fp_line[3])):]

                            #if the variants are the same, writes out to TP vcf file  
                            if fp_sequence == truth_sequence:
                                tp_cor.write(line)
                                tp_read.append([fp_line[1], fp_line[3], fp_line[4]])
                                tp_read.append(left_aligned[0])
                            else:
                                fp_cor.write(line)
                                fp_file_list.append(line) 
    

            else:
                for x in truthset_broken:
                    if int(fp_line[1]) == int(x[0]):
                        if (fp_line[3] == x[1]) & (fp_line[4] == x[2]):
                            tp_cor.write(line)
                            tp_read.append(x)
                        else:
                            fp_cor.write(line)
                            fp_file_list.append(line)            





#opens fn file for reading
with gzip.open(fn_path,"r") as f:
    for line in f:
        if not line.startswith(b'#'):
            fn_line = list(line.decode("utf-8").split('\t'))
            #print(fn_line)

            if {str(fn_line[1])}.issubset(chain.from_iterable(tp_read)):
                for t in tp_read:
                    if (int(t[0]) == int(fn_line[1])):
                        if not (str(t[1]) == str(fn_line[3])) and (str(t[2]) == str(fn_line[4])):
                            fn_cor.write(line) 
            else:
                fn_cor.write(line)

 
            #if not {fn_line[1]}.issubset(chain.from_iterable(vcf_broken)):
                #print("FALSE") 
                #fn_cor.write(line)
                #fn_file_list.append(fn_line)

            #else:
                #print("TRUE")
                #for i, x in enumerate(vcf_broken):
                    #if int(fn_line[1]) == int(x[0]): 
                        #if (fn_line[3] == x[1]) & (fn_line[4] == x[2]): 
                            #if not vcf_file[i] in tp_read:  
                                #tp_cor.write(vcf_file[index])



#close files after writing
fp_cor.close()
tp_cor.close()
fn_cor.close()

