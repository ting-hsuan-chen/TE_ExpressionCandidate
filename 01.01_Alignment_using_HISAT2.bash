#!/bin/bash
#remove adapters by fastq-mcf and check quality of cleaned-up reads
prefix=prefix_of_sample
tRNArRNA=/path/to/tRNA_rRNA.fasta
genome=/path/to/reference_genome.fasta

fastq-mcf -f \
-o ${prefix}_clean_R1.fastq \
-o ${prefix}_clean_R2.fastq \
-l 50 -q 15 -t 0.1 -C 1000000 \
truseq_adapters.fasta \
${prefix}_R1.fastq \
${prefix}_R2.fastq

fastqc ${prefix}_clean_*.fastq -o QC.dir

#build Hisat rRNA indices
hisat2-build ${tRNArRNA} ${tRNArRNA}_hisat2

#build Hisat genome indices
hisat2-build ${genome} ${genome}_hisat2

#align reads to tRNA and rRNA, and then filter for unmapped pairs
hisat2 -x ./${tRNArRNA}_hisat2 \
-1 ${prefix}_clean_R1.fastq \
-2 ${prefix}_clean_R2.fastq \
-S ${prefix}_ncRNA.sam

samtools view -Su -f 12 ${prefix}_ncRNA.sam > ${prefix}_ncRNA.bam

bedtools bamtofastq -i ${prefix}_ncRNA.bam \
-fq ${prefix}_pairs_R1.fastq \
-fq2 ${prefix}_pairs_R2.fastq

#map unmapped reads to Vitis vinifera genome
hisat2 --rna-strandness RF \
--dta -k 100 \
-x Vitis_genome_hisat2 \
-1 ${prefix}_pairs_R1.fastq \
-2 ${prefix}_pairs_R2.fastq \
-S ${prefix}_k100Mm.sam

samtools view -Su ${prefix}_k100Mm.sam |\ samtools sort - ${prefix}_k100Mm
samtools flagstat ${prefix}_k100Mm.bam
samtools index ${prefix}_k100Mm.bam

#sort bam file by name
samtools sort -n -m 4G ${prefix}_k100Mm.bam ${prefix}_k100MmNamesorted.bam
