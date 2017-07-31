#!/bin/bash

# make bwa indices for each reference sequence
for ref in "ribosome_hybrid"
do
    reference="data/reference_sequences/${ref}.fasta"
    bwa index $reference
done

# loop over all samples
for index in {1..105}
do

    echo "Uncompressing index: $index"
    gunzip data/100PE/Index_${index}_*R1*.fastq.gz
    gunzip data/100PE/Index_${index}_*R2*.fastq.gz
        
    echo "Filtering adapter sequences..."
    java -classpath utilities/Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 4 -phred33 data/100PE/Index_${index}_*R1*.fastq data/100PE/Index_${index}_*R2*.fastq ${index}_trimmed_P1.fq ${index}_trimmed_U1.fq ${index}_trimmed_P2.fq ${index}_trimmed_U2.fq ILLUMINACLIP:utilities/Fast-Plast/bin/NEB-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40

    # make reference assemblies of chloroplast, nuclear rDNA, mitochondrial genes for this sample:
    for ref in "ribosome_hybrid"
    do

        reference="data/reference_sequences/${ref}.fasta"
        mkdir data/assemblies/$ref
        output_dir="data/assemblies/${ref}/"
        
        echo "Assembling index: $index"
        
        echo "Mapping reads with bwa..."
        bwa mem $reference -R "@RG\tID:1\tSM:1" ${index}_trimmed_P1.fq ${index}_trimmed_P2.fq > ${index}.sam
        
        echo "Coverting SAM to BAM..."
        samtools view -bS ${index}.sam > ${index}.bam
        
        echo "Sorting reads..."
        samtools sort ${index}.bam ${index}_sorted
        
        echo "Determining read depth..."
        samtools depth ${index}_sorted.bam > ${index}_read_depth.tsv
        
        echo "Making BAM pileup, filtering by phred quality 20..."
        samtools mpileup -Q 20 -Agf $reference ${index}_sorted.bam > ${index}.mpilup

        echo "Generating consensus genotypes..."
        bcftools view -s data/samples.txt -cg ${index}.mpilup > ${index}_temp.vcf
        
        echo "Filtering for read depth >= 10..."
        python scripts/filter_vcf_read_depth.py ${index}_temp.vcf 10

        echo "Generating final FASTA sequence..."
        vcfutils.pl vcf2fq ${index}_temp.vcf.filtered > ${index}.fastq
        seqtk seq -A ${index}.fastq > ${index}.fasta
        
        echo "Cleaning up..."
        rm ${index}_temp.vcf ${index}_temp.vcf.filtered ${index}.mpilup
        rm ${index}_sorted.bam ${index}.bam ${index}.sam ${index}.fastq
        mkdir ${output_dir}${index}
        mv ${index}.fasta ${index}_read_depth.tsv ${output_dir}${index}
    
    done

    echo "Re-compressing index: $index"
    rm ${index}_trimmed_*
    gzip data/100PE/Index_${index}_*R1*.fastq
    gzip data/100PE/Index_${index}_*R2*.fastq

done

# Illumina adapter sequences were removed and the raw sequence reads were quality filtered using Trimmomatic v0.36 (CITE). Reads were trimmed when the average phred quality score in a 10-bp sliding window was less than 20. Reads that were less than 40 bp or that did not survive the filtering process in both forward and reverse directions were excluded.

