# TE_ExpressionCandidate

A computational workflow for the identification of potentially expressed TE loci, which we termed TE expression candidates.


![Image of Workflow](https://github.com/ting-hsuan-chen/TE_ExpressionCandidate/blob/main/Workflow.jpg)
__Figure 1.__ The three-channel workflow for the identification of potentially expressed TE loci



This is a three-channel workflow (Figure 1) mainly comprised of three different tools to quantify TE expression and collect expressed TE loci as follows:

### Sub-pipeline 1
(Corresponding to scripts 01.01, 01.03, and 01.06)

This sub-pipeline only collects TEs obtaining unique-mapping reads that each of these sequencing reads can be traced back to a unique origin in the genome. As illustrated in Figure 1, sequencing reads unmapped to grapevine’s tRNA and rRNA sequences were aligned to 12X PN40024 grapevine reference genome using [HISAT2](https://daehwankimlab.github.io/hisat2/) (Kim et al., 2015a) with the parameters: `<–rna-stradness RF –dtk –k 100>`.  Reads mapping to TEs were than quantified by [htseq-count](https://htseq.readthedocs.io/en/master/count.html) (Anders et al., 2015), in which only uniquely mapped reads were counted. TEs having read count more than 10, which approximates to 5 pairs of read, were collected.

### Sub-pipeline 2
(Corresponding to scripts 01.01, 01.04, and 01.07)

The concept of this sub-pipeline was to collect individual TEs that were aligned with any kind of reads, irrespective of the number of highest quality mapping loci for a given read. Following read alignment using HISAT2 as described in sub-pipeline 1,the [BEDtools](https://bedtools.readthedocs.io/en/latest/) suite (Quinlan and Hall, 2010) was used in read quantification (Figure 1, sub-pipeline 2). The command bedtools coverage generated raw count for TEs while multi-mapping reads matching to n-places (e.g. 10) was recorded n-times (i.e. 10). It also counted number of bases of a TE locus covered by reads (covered bases of TE). Furthermore, this sub-pipeline incorporated bedtools intersect to calculate the number of bases of a read overlapping with an individual TE (mapped bases of read), which were summed for each TE locus to estimate the average read depth of an individual TE’s mapped region (i.e. only the region covered by reads, not the entire annotated feature). This was calculated as follows:
For n reads mapping to a TE locus, and i as integer from 0 to n, f(i) = mapped bases of read i.
Average read depth of an individual TE’s mapped region = <img src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Csum_%7Bi%3D1%7D%5En%20f(i)%20%7D%5Ctext%7BBases%20of%20a%20TE%20locus%20covered%20by%20reads%7D%0A"> .
In order to exclude TE loci that were covered by reads in a sparse and scattered way, a cut-off threshold of the average read depth 5 was adopted in addition to the 10 read count threshold. Examples of the filtering step of this sub-pipeline were illustrated in Figure 2.

(Figure 2)

### Sub-pipeline 3
(Corresponding to scripts 01.02, 01.05, and 01.08)

The third part (Figure 1, sub-pipeline 3) specifically collects TEs had transcription across the boundaries of the element. The software [TEFingerprint](https://github.com/PlantandFoodResearch/TEFingerprint) (Plant and Food Research) was originally designed for identifying unannotated insertions in genomes using paired-end short fragment DNA sequence data. Here it was applied to capture TE loci internally mapped by multi-mapping reads only yet the read mates, as known as danglers, were uniquely aligned to a location near the insertion site (Figure 3). Although TEFingerprint also has the option to use reads sit across the junction of a TE locus, this function was disabled in sub-pipeline 3 as the utilization of htseq-count in sub-pipeline1 has covered this scenario. Sub-pipeline 3 applied the standard TEFingerprint pipeline where reads were mapped against the collection of 223,411 annotated V. vinifera TE sequences using [BWA](https://github.com/lh3/bwa) (Li and Durbin, 2009). Subsequently the mates of TE-mapped reads were aligned to the reference genome (12X PN40024) before calculating read count of dangler clusters. Only clusters containing more than 10 dangler reads were kept to test for the intersection of dangler clusters and annotated TEs using bedtools intersect. The candidates need to show more than 10 dangler reads and more than 10 reads mapping internally (counted by bedtools coverage).

(Figure 3)

### Summarizing TE expression candidates
(Corresponding to the script 01.09)
After excluding TEs did not show enough signal of transcription, potentially expressed TEs from the three sub-pipelines were collected together as a pool of expression candidates.

