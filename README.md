# PROTEINcleaner
a python utility to clean PROTEIN sequences and headers

```
usage: PROTEINcleaner [-h] -i INPUT [-f INPUT_FORMAT] [-o OUTPUT]
                      [-F OUTPUT_FORMAT] [-c Y/y or N/n] [-v Y/y or N/n] [-V]

Cleans invalid bases/residues in PROTEIN sequence file i.e., replaces invalid bases with Xs

Citation: Singh, Abhijeet. PROTEINcleaner: a python utility to clean PROTEIN sequences and headers. ResearchGate 2022, http://dx.doi.org/10.13140/RG.2.2.11468.49283, Available at GitHub: https://github.com/abhijeetsingh1704/PROTEINcleaner

required arguments:
  -i INPUT, --input INPUT
                        input file

options:
  -h, --help            show this help message and exit
  -f INPUT_FORMAT, --input_format INPUT_FORMAT
                        format of input file
                        1 = fasta (default)
                        2 = genbank
  -o OUTPUT, --output OUTPUT
                        output file (default: Clean<time>_<input_file><.ext>)
  -F OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                        format of output fasta file
                        1 = interleaved fasta (default)
                        2 = fasta-2line
                        3 = tab-delimited
  -c Y/y or N/n, --clean_header Y/y or N/n
                        clean and replace special characters in header
                        with underscore "_" (default: N)
  -v Y/y or N/n, --verbose Y/y or N/n
                        print progress to the terminal (default: verbose)
  -V, --version         show program's version number and exit

```

# Citation
Citation: Singh, Abhijeet. PROTEINcleaner: a python utility to clean PROTEIN sequences and headers. ResearchGate 2022, http://dx.doi.org/10.13140/RG.2.2.11468.49283, Available at GitHub: https://github.com/abhijeetsingh1704/PROTEINcleaner
