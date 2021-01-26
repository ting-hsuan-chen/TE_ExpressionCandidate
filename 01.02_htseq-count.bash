#!/bin/bash
#usage: ./01.02_htseq-count.bash -b [name_sorted.bam] -T [TEs.gtf] -G [gene.gtf] -p [prefix_of_sample]
# -b [bam=name_sorted.bam]
# -T [TEgtf=TEs.gtf]
# -G [GENEgtf=gene.gtf]
# -p [prefix=prefix_of_sample]

while getopts b:T:G:p flag
do
    case "${flag}" in
        b) bam=${OPTARG};;
        T) TEgtf=${OPTARG};;
        G) GENEgtf=${OPTARG};;
        p) prefix=${OPTARG};;
    esac
done

#htseq-count --stranded=reverse
#for TEs
htseq-count -f bam \
-t exon \
-i transcript_id \
-s reverse \
-m intersection-nonempty \
$bam \
$TEgtf > ${prefix}_k100Mm_AllSenseTEm_counttable.txt;

#for genes
htseq-count -f bam \
-t exon \
-i transcript_id \
-s reverse \
-m intersection-nonempty \
$bam \
$GENEgtf > ${prefix}_k100Mm_AllSenseGene_counttable.txt;
