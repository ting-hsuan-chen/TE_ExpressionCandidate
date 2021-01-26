#!/bin/bash
#usage: ./01.06_bwa.bash -g [reference_genome.fasta] -f [TE.fasta] -p [prefix_of_sample]
# -g [genome=reference_genome.fasta]
# -f [TEfasta=TE.fasta]
# -p [prefix=prefix_of_sample]

while getopts g:f:p: flag
do
    case "${flag}" in
        g) genome=${OPTARG};;
        f) TEfasta=${OPTARG};;
        p) prefix=${OPTARG};;
    esac
done

#Build index
bwa index -a bwtsw $genome;
bwa index -a is $TEfasta;

#Align read pairs against sequences of annotated TEs
bwa mem -t 8 \
    AllTEsExpanded_curated.fa \
    ${prefix}_pairs_R1.fastq \
    ${prefix}_pairs_R2.fastq \
    | samtools view -Su - \
    | samtools sort - -o ${prefix}_TEmapped.bam;


#index the resulting bam file
samtools index /${prefix}_TEmapped.bam;
