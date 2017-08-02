#!/bin/bash

reference="reference_genomes/KX118606.gb"
input_dir="assemblies/bwa/chloroplast/"
output_dir="annotated_plastomes/"

for index in 33 34 #35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58
do
    # plann.pl -reference gbfile.gb -fasta denovoplastome.fasta -out outfile [-organism "Genus species"] [-sample samplename]
    perl utilities/plann/plann.pl -reference $reference -fasta ${input_dir}${index}/${index}.fasta -out ${output_dir}${index}
    
    # make genbank files
    # tbl2asn -t template.sbt -i 33.fsa -V vb
done
