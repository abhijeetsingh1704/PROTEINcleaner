#!/usr/bin/env python3
###################################
# File      - PROTEINcleaner
# Modified  - sat 28 maj 2022 15:37:18 CEST
# Sign      - Abhijeet
###################################
program = 'PROTEINcleaner'
version = ': version (1.0)\nAuthor: Abhijeet Singh <abhijeetsingh.aau@gmail.com>'
citation = '''\n\nCitation: Singh, Abhijeet. PROTEINcleaner: a python utility to clean PROTEIN sequences and headers. ResearchGate 2022, http://dx.doi.org/10.13140/RG.2.2.11468.49283, Available at GitHub: https://github.com/abhijeetsingh1704/PROTEINcleaner'''
###################################
import sys
import os
import datetime
import subprocess
from collections import defaultdict
import argparse
import string
import re
#
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Biopython missing, attempting to install...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython>=1.79"])
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
#
###################################
# parser
parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description='Cleans invalid bases/residues in PROTEIN sequence file i.e., replaces invalid bases with Xs' + citation)
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--input", dest='input', required=True, type=str, help='input file')
optional.add_argument("-f", "--input_format", dest='input_format', required=False, type=int, default=1, help='format of input file\n1 = fasta (default)\n2 = genbank\n')
optional.add_argument("-o", "--output", dest='output', required=False, type=str, help='output file (default: Clean<time>_<input_file><.ext>)')
optional.add_argument("-F", "--output_format", dest='output_format', required=False, type=int, default=1, help='format of output fasta file\n1 = interleaved fasta (default)\n2 = fasta-2line\n3 = tab-delimited\n')
optional.add_argument("-c", "--clean_header", dest='clean_header', required=False, metavar='Y/y or N/n', default='N', type=str, help='clean and replace special characters in header\nwith underscore "_" (default: N)')
optional.add_argument("-v", "--verbose", dest='verbose', required=False, metavar='Y/y or N/n', default='Y', type=str, help='print progress to the terminal (default: verbose)')
optional.add_argument('-V', '--version', action='version', version='%(prog)s'+ str(version))
parser._action_groups.append(optional)
args = parser.parse_args()
###################################
#
def filenames(filename):
    absname = os.path.abspath(filename)
    basename = os.path.basename(filename)
    file_name = os.path.splitext(basename)[0]
    file_ext = os.path.splitext(basename)[1]
    return absname, basename, file_name, file_ext
###################################
#
date_var_1 = datetime.datetime.now()
date_var_2 = date_var_1.strftime("%H%M%S")
###################################
# 
input_args = args.input
output_args = args.output
input_format_args = args.input_format
output_format_args = args.output_format
verbosity_args = args.verbose.upper()
cleaner_args = args.clean_header.upper()
#
if cleaner_args == "Y":
    clean_word = "YES"
elif cleaner_args == "N":
    clean_word = "NO"
###################################
#
input_file = filenames(input_args)[0]
###################################
# optional options input format
if type(input_format_args) == int:
    if input_format_args == 1:
        input_format = "fasta"
    if input_format_args == 2:
        input_format = "genbank"
    if input_format_args >= 3:
        sys.exit("[ERROR]: input format flag requires interger value\n\t use -f or --input_format with value 1 or 2")
else:
    sys.exit("[ERROR]: input format flag requires interger value\n\t use -f or --input_format with value 1 or 2")
###################################
# optional options output format
if type(output_format_args) == int:
    if output_format_args == 1:
        output_format = "fasta"
        output_ext = ".fasta"
    if output_format_args == 2:
        output_format = "fasta-2line"
        output_ext = ".fasta"
    if output_format_args == 3:
        output_format = "tab"
        output_ext = ".tab"
    if output_format_args > 3:
        sys.exit("[ERROR]: output format flag requires interger value\n\t use -F or --output_format with value 1, 2 or 3")
else:
    sys.exit("[ERROR]: output format flag requires interger value\n\t use -F or --output_format with value 1, 2 or 3")
###################################
# output
if output_args is not None:
    output_var = os.path.abspath(output_args)
    output_file = output_var + output_ext
    output_obj = open(output_file, "w")
else:
    output_var = os.path.join(os.getcwd(), "Clean" + date_var_2 + "_" + filenames(input_file)[2])
    output_file = output_var + output_ext
    output_obj = open(output_file, "w")
###################################
# change objects
Protein_AA = list("ACDEFGHIKLMNPQRSTVWYX")
numbers = str("0123456789")
All_Characters = list(string.ascii_uppercase + string.punctuation + numbers) #.replace(".","")
Remove_Characters = [x for x in All_Characters if x not in Protein_AA]
# assign start values
input_num = 0
aa_change = 0
aa_change_seq_number = 0
aa_change_seq_type = 0
aa_change_all = 0
invalid_list = []
###################################
# parse input file
query = SeqIO.index(input_file, input_format)
for key, values in query.items():
    PROTEIN_id = values.id
    PROTEIN_description = values.description.replace(values.id, "").strip(" ?")
    PROTEIN_sequence = str(values.seq.upper())
    #
    if cleaner_args == "Y":
        PROTEIN_description_2 = str(PROTEIN_description.replace(" ", "_"))
        for BadChr in string.punctuation:
            PROTEIN_description_2 = PROTEIN_description_2.replace(BadChr, "_")
            PROTEIN_description_2 = re.compile(r'_{2,}').sub('_', PROTEIN_description_2)
    #
    for AA in Remove_Characters:
        if AA in PROTEIN_sequence:
            invalid_list.append(AA)
            aa_change_seq_number = str(PROTEIN_sequence).count(AA)
            PROTEIN_sequence = PROTEIN_sequence.replace(AA, "N")
            aa_change_seq_type += int(aa_change_seq_number)
    #
    if verbosity_args == "Y":
        if cleaner_args == "Y":
            if output_format_args == 3:
                print(PROTEIN_id + "_" + PROTEIN_description_2 + "\t" + PROTEIN_sequence)
            elif output_format_args == 1 or 2:
                print(">" + PROTEIN_id + "_" + PROTEIN_description_2 + "\n" + PROTEIN_sequence)
        #        
        elif cleaner_args == "N":
            if output_format_args == 3:
                print(PROTEIN_id + " " + PROTEIN_description + "\t" + PROTEIN_sequence)
            elif output_format_args == 1 or 2:
                print(">" + PROTEIN_id + " " + PROTEIN_description + "\n" + PROTEIN_sequence)
    # cleaned sequences
    if cleaner_args == "Y":
        PROTEIN_header = PROTEIN_id + "_" + PROTEIN_description_2
        final_seq = (SeqRecord(Seq(str(PROTEIN_sequence)), id=PROTEIN_header, description="", name=""))
    else:
        if output_format_args == 3:
            PROTEIN_header = PROTEIN_id + " " + PROTEIN_description # 
            final_seq = (SeqRecord(Seq(str(PROTEIN_sequence)), id=PROTEIN_header, description="", name=""))
        elif output_format_args == 1 or 2:
            final_seq = (SeqRecord(Seq(str(PROTEIN_sequence)), id=PROTEIN_id, description=PROTEIN_description, name=""))
    # write output
    SeqIO.write(final_seq, output_obj, output_format)
    #
    input_num += 1
#    
aa_change_all += int(aa_change_seq_type)
#
output_obj.close()
###################################
#
if not invalid_list:
    invalid_list = "none"
else:
    invalid_list = str(",".join(list(set(invalid_list))))
###################################
# Print info
print("#" * 60)
print("[Program]\t: " + program)
print("[Date]\t\t: "+ date_var_1.strftime("%Y-%m-%d %H:%M:%S"))
print("[Valid AAs]\t: " + str(",".join(Protein_AA)))
print("[Input file]\t: "+ input_file)
print("\t\t  |_[Invalid AAs]: " + invalid_list)
print("[Output file]\t: "+ output_file)
print("\t\t  |_[Output Seqs]: " + str(aa_change_all) + " AA change(s) in " +  str(input_num) + " sequence(s)")
print("[Header cleaned]: "+ clean_word)
print("#" * 60)
##################################
# end of script