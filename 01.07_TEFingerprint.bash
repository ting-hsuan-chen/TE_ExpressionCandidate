#!/bin/bash
#tef preprocess (no clipped reads used)
mydir=/home/user
source activate tefingerprint
export PATH=$HOME/tools/TEFingerprint-master/applications/:$PATH
cd /home/user/tools/TEFingerprint-master/

for i in {01..39}; do
tef preprocess \
$mydir/${i}_TEmapped.bam \
--reference $mydir/Vitis_genome.fasta \
--output $mydir/preprocess/${i}_danglers.bam \
--exclude-tails \
--tempdir $mydir/temp \
--threads 8;
done

#generate unique TE family lists
awk -vOFS='\t' '{print $1}' $mydir/TEfamilies_named.lst \
| sort | uniq > $mydir/TEfamilies.txt

#tef compare
#threshold for TE-related dangler > 10 reads, so the minimum value is set to 11 (-m 11)
tef compare $mydir/preprocess/01_danglers.bam $mydir/preprocess/02_danglers.bam \
$mydir/preprocess/03_danglers.bam $mydir/preprocess/04_danglers.bam \
$mydir/preprocess/05_danglers.bam $mydir/preprocess/06_danglers.bam \
$mydir/preprocess/07_danglers.bam $mydir/preprocess/08_danglers.bam \
$mydir/preprocess/09_danglers.bam $mydir/preprocess/10_danglers.bam \
$mydir/preprocess/11_danglers.bam $mydir/preprocess/12_danglers.bam \
$mydir/preprocess/13_danglers.bam $mydir/preprocess/14_danglers.bam \
$mydir/preprocess/15_danglers.bam $mydir/preprocess/16_danglers.bam \
$mydir/preprocess/17_danglers.bam $mydir/preprocess/18_danglers.bam \
$mydir/preprocess/19_danglers.bam $mydir/preprocess/20_danglers.bam \
$mydir/preprocess/21_danglers.bam $mydir/preprocess/22_danglers.bam \
$mydir/preprocess/23_danglers.bam $mydir/preprocess/24_danglers.bam \
$mydir/preprocess/25_danglers.bam $mydir/preprocess/26_danglers.bam \
$mydir/preprocess/27_danglers.bam $mydir/preprocess/28_danglers.bam \
$mydir/preprocess/29_danglers.bam $mydir/preprocess/30_danglers.bam \
$mydir/preprocess/31_danglers.bam $mydir/preprocess/32_danglers.bam \
$mydir/preprocess/33_danglers.bam $mydir/preprocess/34_danglers.bam \
$mydir/preprocess/35_danglers.bam $mydir/preprocess/36_danglers.bam \
$mydir/preprocess/37_danglers.bam $mydir/preprocess/38_danglers.bam \
$mydir/preprocess/39_danglers.bam \
   	-f $(< $mydir/TEfamilies.txt) \
   	-m 11 \
   	-e 600 \
   	-q 30 \
   	-t 8 \
   	> $mydir/compare/AllRep_comparison.gff
#use "Excel" to open the file and define tab, comma, semicolon and "=" as separators, #and then save the file as tab-separated text file.

#capture TE-related danglers overlapping with annotated TEs using “bedtools intersect”
#!/bin/bash
cd ~
bedtools intersect -a ./compare/AllRep_comparison.gff \
-b AllTEsExpanded_curated_V3.gtf \
-wa \
-wb \
> AllRep_comp_intersect.gtf
#use "Excel" to open the output file and define tab, comma, semicolon, space and "=" as #separators, and then save it as tab-separated text file: AllRep_comp_intersect_V4.gtf.
