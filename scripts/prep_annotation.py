# code used to preprocess the bins and create the directory structure 
# for output

# make working directory stucture
#  -Annotation_Output
#    -Gene_calls
#    -Protein_ortho_out
#    -Annotation_results

import os
import sys
from Bio import SeqIO
import glob

def length_generator(bin_file):
    for record in SeqIO.parse(bin_file,"fasta"):
        yield(len(record.seq))

def check_len(total_len, min_len):
    if total_len > min_len:
        return(True)
    else:
        return(False)

path_to_bins = sys.argv[1].strip("/")
MIN_LEN = 20000

# make directories
os.makedirs("Passing_Bins")
os.makedirs("Gene_Calls")
os.makedirs("Protein_Ortho_Out")
os.makedirs("Annotation_Results")

passing_bins_dir = "Passing_Bins/"

# filter bins by size
bins = [x for x in glob.glob(path_to_bins + '/*')]

for bin_file in bins:
    temp_gen = length_generator(bin_file)
    total_len = 0
    for length in temp_gen:
        total_len += length

    if check_len(total_len,MIN_LEN):
        print("pass")
        os.symlink(os.getcwd()+"/"+bin_file,passing_bins_dir + bin_file.split("/")[-1])
    else:
        print("fail")    