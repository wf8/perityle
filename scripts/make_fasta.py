#! /usr/bin/python
"""
Combine individual fasta files for each specimen into a single fasta file
with the correct sequence id (Taxon_name_CHY00).
"""

import csv
import os
import sys
from os import listdir
from os.path import isfile, join

index_file = "perityle_specimen_indices.csv"
input_dir = "data/assemblies/" 
output_dir = "data/unaligned/"
genes = ["ATP1", "ATP4", "ATP6", "ATP8", "ATP9", "chloroplast", "COB", "COXI", "COXIII", "mitochondrion", "NAD3", "NAD4L", "NAD5", "NAD6", "NAD9", "ribosome", "RPL16", "RPS12", "RPS13", "RPS4", "RRN18", "RRN26", "RRN5", "ribosome_hybrid"]

indices = {}

with open(index_file, 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    for row in csvreader:
        if row[1] != '' and row[1] != 'sequence id':
            indices[ row[1] ] = row[5] + ' ' + row[6] + ' ' + row[0]
            indices[ row[1] ] = indices[ row[1] ].replace(' ', '_').replace('.', '')

for gene in genes:
    with open(output_dir + gene + '.fasta', 'w') as fout:
        for key in indices:
            fout.write('>' + indices[ key ] + '\n')
            # data/assemblies/ATP1/1/1.fasta
            with open(input_dir + gene + '/' + key + "/" + key + ".fasta", 'r') as seqfile:
                seq = seqfile.readlines()[1]
                fout.write(seq)
