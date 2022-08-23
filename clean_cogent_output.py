#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""concatenates Cogent output with unassigned input to bring together Cogent's 'reconstructed transcriptome'"""

from __future__ import print_function

import os
import sys
import argparse
from Bio import SeqIO


def main(fasta, cogent_output_dir, fileout):
    # get a list of sequences that were partitioned
    partitioned_seq_ids = []
    partitions = os.listdir(cogent_output_dir)
    for partition in partitions:
        # open each fasta file and record the sequence IDs
        for seq in SeqIO.parse('{}/{}/in.fa'.format(cogent_output_dir, partition), 'fasta'):
            partitioned_seq_ids.append(seq.id)

    # get all the hq sequences
    not_partitioned = []
    for seq in SeqIO.parse(fasta, 'fasta'):
        # subtract partitioned from all to get not-partitioned
        if seq.id not in partitioned_seq_ids:
            not_partitioned.append(seq)

    # check output, and log any trouble / add all sequences that ran correctly
    all_worked = True
    partitioned = []
    for partition in partitions:
        expected_result = '{}/{}/cogent2.renamed.fasta'.format(cogent_output_dir, partition)
        try:
            for seq in SeqIO.parse(expected_result, 'fasta'):
                partitioned.append(seq)
        except IOError:
            print("Warning: {} not found".format(expected_result), file=sys.stderr)
            all_worked = False
    if not all_worked:
        print("""For cases where the output file is missing, the cogent Readme suggests rerunning 
'reconstruct_contig.py' with, e.g., --nx_cycle_detection -k 40""")
    # output the original non partitioned sequences, and where they were before cogent, and the further collapsed
    # cogent output for the partitioned sequences
    out = not_partitioned + partitioned
    with open(fileout, 'w') as f:
        SeqIO.write(out, f, 'fasta')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='reorganize output of cogent to produce one fasta for the '
                                                 'whole transcriptome')
    parser.add_argument('-f', '--fasta', required=True, type=str, help='fasta file that cogent received as input')
    parser.add_argument('-c', '--cogent_output_dir', required=True, type=str, help='output directory for cogent')
    parser.add_argument('-o', '--out', required=True, type=str, help='output file name')

    args = parser.parse_args()
    main(args.fasta, args.cogent_output_dir, args.out)


