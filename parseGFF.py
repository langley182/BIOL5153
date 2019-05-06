#! /usr/bin/env python3

# Create a GFF parser that takes as input the watermelon genome and the gff report,
# and outputs a file that contains all the separate genes


# 4/3/19 - New instructions - add the argparse stuff, and have 2 positional arguments - the gff and the fasta files. 

# specify the input files
gff_file = 'watermelon.gff'
fasta_file = 'watermelon.fsa'

# open the fasta file
fasta = open(fasta_file, 'r')

# open the GFF file
gff = open(gff_file, 'r')

# read the GFF file, line by line
for line in gff:
    # Skip blank lines
    if not line.strip()
    
    # remove line breaks
    line = line.rstrip('\n')
    # print (line) 

    # split each line on the tab character
    sequence, source, feature, begin, end, length, strand, phase, attributes = line.split('\t')
    
    # print the DNA sequence for this feature

    # extract the DNA sequence from the genome for this feature

    # Calculate the GC content 
    

# end

gff.close()