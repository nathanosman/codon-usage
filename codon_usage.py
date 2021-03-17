"""
This script takes as an input a FASTA file that contains the nucleotide sequence and calculates the codon usage for each sequence
This is a Python port of the script (originally written in Java, and later in VBA) used for this article: "Palanisamy N, Osman N, Ohnona F, et al. Does antiretroviral treatment change HIV-1 codon usage patterns in its genes: a preliminary bioinformatics study. AIDS Res Ther. 2017;14(1):2. Published 2017 Jan 7. doi:10.1186/s12981-016-0130-y"

Author:
    Nathan Osman (nathan.osman@mail.mcgill.ca)

Version:
    v2.0 (2021-03-17)

Example:
    python codon_usage.py sequence_file.fasta


Ideas for improvements:
    - Add csv option
    - Add svg option
    - Add support for other file formats
    - Add support for different genetic codes

"""

import argparse
import pandas as pd
import sys
from Bio import SeqIO
from Bio.Data import CodonTable as GC

# Set the parser that will take care of the script's arguments
parser = argparse.ArgumentParser()
parser.add_argument('infile', help='Path to the input file', nargs='?')
# parser.add_argument('outfile', help='Path to the output file', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
args = parser.parse_args()

# Gather the sequences in a dictionnary of SeqRecords
record_dict = SeqIO.index(args.infile, 'fasta')

# Create a pandas MultiIndex containing the different codons for each amino-acid, using the standard genetic code from the Bio.Data package
gc_dict = GC.unambiguous_dna_by_id[1].forward_table
gc_index = pd.MultiIndex.from_tuples([(codon, aa) for codon, aa in gc_dict.items()], names=['aa', 'codon'])

# Create a pandas DataFrame to contain all the frequencies from all the sequences
df = pd.DataFrame(index=gc_index)

# For each sequence in the dictionnary, create a pandas Series containing the codon usage and append to the dataframe
for record in record_dict.values():
    s = pd.Series(index=gc_index.get_level_values('codon'), name=record.id, dtype=int)
    # Go through the sequence codon by codon and increment the number of the corresponding codon in the series
    for i in range(len(record)//3):
        codon = str(record.seq[3*i:3*i+3])
        s.loc[codon] += 1
    
    # Join the series to the dataframe
    df = df.join(s)

print(df)
