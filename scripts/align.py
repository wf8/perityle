#! /usr/bin/python

import csv
from Bio import Entrez
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from time import sleep

genes = ["ATP1", "ATP4", "ATP6", "ATP8", "ATP9", "chloroplast", "COB", "COXI", "COXIII", "mitochondrion", "NAD3", "NAD4L", "NAD5", "NAD6", "NAD9", "ribosome", "RPL16", "RPS12", "RPS13", "RPS4", "RRN18", "RRN26", "RRN5"]

for i, gene in enumerate(genes):
    print("Aligning " + gene + " with MAFFT...")
    mafft_cline = MafftCommandline(input="data/unaligned/" + genes[i] + ".fasta")
    mafft_cline.set_parameter("--auto", True)
    print(str(mafft_cline))
    stdout, stderr = mafft_cline()

    print("Writing " + gene + " alignment to FASTA file...")
    with open("data/aligned/" + genes[i] + ".fasta", "w") as handle:
        handle.write(stdout)

print("Done.")
