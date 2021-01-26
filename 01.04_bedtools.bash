#!/bin/bash
#usage: ./01.04_bedtools.bash -T [TEs.gtf] -p [prefix_of_sample]
# -T [TEgtf=TEs.gtf]
# -p [prefix=prefix_of_sample]

while getopts b:T:G:p flag
do
    case "${flag}" in
        T) TEgtf=${OPTARG};;
        p) prefix=${OPTARG};;
    esac
done

#capture F1R2 read pairs
bamtools filter -script filter_F1R2.js -in ${prefix}_k100Mm.bam -out ${prefix}_k100Mm_F1R2.bam;
samtools index ${prefix}_k100Mm_F1R2.bam;

#capture F2R1 read pairs
bamtools filter -script filter_F2R1.js -in ${prefix}_k100Mm.bam -out ${prefix}_k100Mm_F2R1.bam;
samtools index ${prefix}_k100Mm_F2R1.bam;

#count reads mapping to TEs in sense direction using “bedtools coverage”:
#for TEs in forward orientation, use F2R1 read pairs;
#for TEs in reverse orientation, use F1R2 read pairs.
bedtools bamtobed -split -i ${prefix}_k100Mm_F2R1.bam > ${prefix}_k100Mm_F2R1_split.bed;
bedtools bamtobed -split -i ${prefix}_k100Mm_F1R2.bam > ${prefix}_k100Mm_F1R2_split.bed;
bedtools coverage -a $TEgtf \
-b ${prefix}_k100Mm_F2R1_split.bed \
> ${prefix}_BedCov_forwardTE_sense.txt
bedtools coverage -a AllTEsExpanded_curated_Reverse_V2.gtf \
-b ${prefix}_k100Mm_F1R2_split.bed \
> ${prefix}_BedCov_reversedTE_sense.txt;

cat ${prefix}_BedCov_forwardTE_sense.txt ${prefix}_BedCov_reversedTE_sense.txt \
| sort -k1 \
> ${prefix}_BedCov_TE_sense.txt;


#average read depth = ( Sum of MapReadBase )/( MapTEbase )
#Sum of MapReadBase = sum of mapped bases of the TE-mapping reads
#MapTEbase = sum of the mapped bases of the corresponding TEs

#calculate MapReadBase of each read overlapping with each TE using “bedtools intersect”
bedtools intersect -a $TEgtf \
-b ${prefix}_k100Mm_F2R1_split.bed \
-wo > ${prefix}_BedIntersect_forwardTE_sense.txt;
bedtools intersect -a $TEgtf \
-b ${prefix}_k100Mm_F1R2_split.bed \
-wo > ${prefix}_BedIntersect_reversedTE_sense.txt;

cat ${prefix}_BedIntersect_forwardTE_sense.txt ${prefix}_BedIntersect_reversedTE_sense.txt \
| sort -k1 >  ${prefix}_BedIntersect_TE_sense.txt;


#establish the lists of TEs having sense reads and TEs not having sense reads
awk '{gsub("\"", "");print}' $TEgtf > temp.gtf

awk -F"\t" '{print $9}'  ${prefix}_BedIntersect_TE_sense.txt \
| uniq > ${prefix}_TEhaveSenseRead.txt
awk -F'\t' 'NR==FNR {c[$1]++;next}; c[$9]==0'  ${prefix}_TEhaveSenseRead.txt temp.gtf \
| awk -F"\t" '{print}' > ${prefix}_TEhaveNoSenseRead.txt;
done

rm -f temp.gtf

#calculate “Sum of MapReadBase” for each TEs having sense reads
#required information in the file “${i}_BedIntersect_TE_sense.txt”:
#MapReadBase (number of mapped bases of each TE-mapping read): the last column
touch ${prefix}_temp01.txt;

cat ${prefix}_TEhaveSenseRead.txt | while read line; do
grep "$line" ${prefix}_BedIntersect_TE_sense.txt > ${prefix}_temp02.txt
awk '{sum+=$NF}END{print sum}' ${prefix}_temp02.txt >> ${prefix}_temp01.txt; done;

# calculate “average read depth”
#required information:
#(1) TE id in the file “${i}_BedCov_TE_sense.txt”: 9th column
#(2) TE id in the file “${i}_OverlapBasePair_Read.txt”: 1st column
#(3) Sum of MapReadBase: will be in the 14th column of “${i}temp03.txt”
#(4) MapTEbase: will be in the 11th column of “${i}temp03.txt”

paste -d "\t" ${prefix}_TEhaveSenseRead.txt ${prefix}_temp01.txt > ${prefix}_OverlapBasePair_Read.txt;

join -1 9 -2 1 -t $'\t' \
<(sort -k9 ${prefix}_BedCov_TE_sense.txt) \
<(sort -k1 ${prefix}_OverlapBasePair_Read.txt) \
| awk -F'\t' '{OFS="\t"; print $2,$3,$4,$5,$6,$7,$8,$9,$1,$10,$11,$12,$13,$14}' - \
> ${prefix}_temp03.txt;

awk -F'\t' '{depth=$14/$11; print depth}' ${prefix}_temp03.txt \
| paste ${prefix}_temp03.txt - \
> ${prefix}_BedCov_TEhaveSenseRead.txt;

awk -F"\t" '{OFS="\t"; TElength=$5-$4+1; print $1,$2,$3,$4,$5,$6,$7,$8,$9,0,0,
TElength,0,0,0}' ${prefix}_TEhaveNoSenseRead.txt \
> ${prefix}_BedCov_TEhaveNoSenseRead.txt;

cat ${prefix}_BedCov_TEhaveSenseRead.txt ${prefix}_BedCov_TEhaveNoSenseRead.txt \
| sort -k1 \
> ${prefix}_BedCov_OverlapBasePair_senseRead.txt;

rm -f ${prefix}_temp*.txt
rm -f ${prefix}_TEhave*.bed;
