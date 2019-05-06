#! /usr/bin/env python3

import argparse
import csv
import re
from Bio import SeqIO
from collections import defaultdict

def get_args():
    # create an argument parser object
    parser = argparse.ArgumentParser(description = 'Parses a GFF file')
    
    # add positional argument for the input position in the Fibonacci sequence
    parser.add_argument("gff_file", help="GFF3 formated file")
    parser.add_argument("fasta_file", help="FASTA file corresponding to the GFF3 file")

    # parse the arguments
    return parser.parse_args()

def parse_fasta():
    genome_object = SeqIO.read(args.fasta_file, 'fasta')
    return genome_object.seq
    
def parse_gff(genome):
   
   # Dictionary to hold CDS sequences, key = gene name, values = dict #2, where key = exon number
   # value = the sequence
    coding_seqs = defaultdict(dict)
   
    # Open, read, and parse GFF file
    with open(args.gff_file,'r') as gff:

        # Create a csv reader object
        reader = csv.reader(gff, delimiter = '\t')
        
        for line in reader:
            # skip blak lines
            if not line:
                continue
            
            # Skip commented lines
            #elif re.match('^#', line):
             #   continue

            else: 
                begin = int(line[3])-1
                end   = int(line[4])
                strand = line[6]
                feature_type = line[2]
                attributes = line[8]
                species = line[0]

                feature_sequence = genome[begin:end]

                if(strand == '-'):
                    feature_sequence = rev_comp(feature_sequence)

                # Calculate GC content of feature
                gc_content = gc(feature_sequence)

                if(feature_type == 'CDS'):
                    # Split attributes field into its separate parts, to get the gene inro
                    exon_info = attributes.split( ' ; ')

                    # Extract the exon number
                    gene_name = exon_info[0].split()[1]

                    # Exctract the gene name if it is an exon
                    #if(exon_info[0].split()[2]):
                    if len(exon_info[0].split()) > 2:
                        # Exctact the exon number
                        exon_number = exon_info[0].split()[-1]
                        
                        if gene_name in coding_seqs:
                            # Store the coding seq for this exon
                            coding_seqs[gene_name][exon_number] = feature_sequence

                        else:
                            # First time encountering this gene; declare the dict for seq
                            coding_seqs[gene_name] = {}
                            # Store the coding seq for this exon
                            coding_seqs[gene_name][exon_number] = feature_sequence

                    #else:
                        #continue
                        # Print the sequence in FASTA format
                        #print('>' + species.replace(' ','_') + '_' + gene_name)
                        #print(feature_sequence)

    # Done reading the GFF file, loop over coding_seqs to  print the CDS sequence
    # gene = gene name, exons = dict of exon sequences (key = exon num, value = exon seq)
    for gene, exons in coding_seqs.items():
        # Print the fasta header for each gene
        print('>' + species.replace(' ','_') + '_' + gene)

        # Make a variable that will hold the concatenated  CDS sequence
        cds_for_for_this_gene = ''

        # Loop over all exons in the gene
        # Needs to be sorted!
        for exon_num, exon_seq in sorted(exons.items()):
            cds_for_for_this_gene += exon_seq
        
        # Print the CDS for this gene
        print(cds_for_for_this_gene)

def rev_comp(seq):
    return seq.reverse_complement()

def gc(seq):
    seq = seq.upper()
    count_of_G = seq.count('G')
    count_of_C = seq.count('C')

    return (count_of_G + count_of_C) / len(seq)

def main():
    genome = parse_fasta()
    parse_gff(genome)
    
# Get the arguments before calling main
args = get_args()


# Execute the program by calling main
if __name__ == "__main__":
	main()
