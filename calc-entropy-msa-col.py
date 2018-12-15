# This script will calculate Shannon entropy from a MSA

import os
import sys
import warnings
import traceback


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Compute per residue Shannon entropy of a Multiple Sequence Alignment.')

        parser.add_argument('-a',
                            '--alignment',
                            action='store',
                            required=True,
                            help='The multiple sequence alignment (MSA) for protein 1 in any of the formats supported by Biopython\'s AlignIO.')
        #parser.add_argument('-a2',
        #                    '--alignment',
        #                    action='store',
        #                    required=True,
        #                    help='The multiple sequence alignment (MSA) for protein 2 in any of the formats supported by Biopython\'s AlignIO.')
        parser.add_argument('-f',
                            '--alnformat',
                            action='store',
                            default='fasta',
                            help='Specify the format of the input MSA to be passed in to AlignIO.')
        parser.add_argument('-v',
                            '--verbose',
                            action='count',
                            default=0,
                            help='Verbose behaviour, printing parameters of the script.')
        parser.add_argument('-m',
                            '--runningmean',
                            action='store',
                            type=int,
                            default=0,
                            help='Return the running mean (a.k.a moving average) of the MSAs Shannon Entropy. Makes for slightly smoother plots. Providing the number of points to average over switches this on.')
        parser.add_argument('--makeplot',
                            action='store_true',
                            help='Plot the results via Matplotlib.')
    except:
        print ("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()



def parseMSA(msa, alnformat, verbose):
    """Parse in the MSA file using Biopython's AlignIO"""

    from Bio import AlignIO
    alignment = AlignIO.read(msa, alnformat)

    #---check the input such that each alignment has the same number of residues
    seq_lengths_list = []
    for record in alignment:
       seq_lengths_list.append(len(record))

    seq_lengths = set(seq_lengths_list)

    if verbose > 0: print("Alignment length is:" + str(list(seq_lengths)))

    if len(seq_lengths) != 1:
        sys.stderr.write("Your alignment lengths aren't equal. Check your alignment file.")
        sys.exit(1)

    index = range(1, list(seq_lengths)[0]+1)

    return alignment, list(seq_lengths), index

##################################################################
# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i
# Gaps and N's are included in the calculation
##################################################################

def shannon_entropy(list_input):
    """Calculate Shannon's Entropy per column of the alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)"""

    import math
    unique_seq = set(list_input)
    print(unique_seq)
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for seq in unique_seq:
        n_i = list_input.count(seq) # Number of residues of type i
        P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)

    sh_entropy = -(sum(entropy_list))

    return sh_entropy


def shannon_entropy_list_msa(alignment):
    """Calculate Shannon Entropy across the whole MSA"""

    shannon_entropy_list = []  #initialize empty list: shannon entropy for a column
    
    for col_no in range(len(list(alignment[0]))):
        list_input = list(alignment[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))

    return shannon_entropy_list


def mutual_entropy(list_input1, list_input2):

    import math
    unique_seq = set(list_input1, list_input2)



def mutual_entropy_list_msa(alignment):
    """Calculate mutual entropy """

    mutual_entropy_list = []

    

    

def plot(index, sel, verbose):
    """plot via matplotlib to visualize"""
    import matplotlib.pyplot as plt

    if verbose > 0: print("Plotting data...")

    plt.plot(index, sel)
    plt.xlabel('MSA Position Index', fontsize=12)
    plt.ylabel('Shannon Entropy', fontsize=12)

    plt.show()


def main():
    """Compute Shannon Entropy from a provided MSA."""

    # Parse arguments
    args = parseArgs()

    # Convert object elements to standard variables for functions
    msa = args.alignment
    alnformat = args.alnformat
    verbose = args.verbose
    makeplot = args.makeplot

    #----call functions to print output

    alignment, seq_lengths, index = parseMSA(msa, alnformat, verbose)
    #sel = shannon_entropy_list_msa(alignment), shannon_entropy(alignment)
    sel = shannon_entropy_list_msa(alignment)

    #if runningmean > 0:
    #    sel = running_mean(sel, runningmean)

    if makeplot:
        plot(index, sel, verbose)

    if verbose > 0:
        print("Index" + '\t' + "Entropy")
        
    for c1, c2 in zip(index, sel):
        print(str(c1) + '\t' + str(c2))


if __name__ == '__main__':
    main()
