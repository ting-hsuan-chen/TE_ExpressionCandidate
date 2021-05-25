# TE_ExpressionCandidate
### A computational workflow for the identification of potentially expressed TE loci, which we termed TE expression candidates.
This is a three-channel workflow (Figure 1) mainly comprised of three different tools to quantify TE expression and collect expressed TE loci as follows:
<br>

<img src="https://github.com/ting-hsuan-chen/TE_ExpressionCandidate/blob/main/Workflow.jpg" width="800">

__Figure 1.__ The three-channel workflow for the identification of potentially expressed TE loci
<br>
<br>

### Sub-pipeline 1
(Corresponding to scripts 01.01_HISAT2.bash, 01.03_htseq-count.bash, and 01.06_filter_subpipeline_1.r)

This sub-pipeline only collects TEs obtaining unique-mapping reads that each of these sequencing reads can be traced back to a unique origin in the genome. As illustrated in Figure 1, sequencing reads unmapped to tRNA and rRNA sequences were aligned to the reference genome using [HISAT2](https://daehwankimlab.github.io/hisat2/) (Kim et al., 2015a) with the parameters: `<–rna-stradness RF –dtk –k 100>`.  Reads mapping to TEs were than quantified by [htseq-count](https://htseq.readthedocs.io/en/master/count.html) (Anders et al., 2015), in which only uniquely mapped reads were counted. TEs having read count more than 10, which approximates to 5 pairs of read, were collected.

### Sub-pipeline 2
(Corresponding to scripts 01.01_HISAT2.bash, 01.04_bedtools.bash, and 01.07_filter_subpipeline_2.r)

The concept of this sub-pipeline was to collect individual TEs that were aligned with any kind of reads, irrespective of the number of highest quality mapping loci for a given read. Following read alignment using HISAT2 as described in sub-pipeline 1,the [BEDtools](https://bedtools.readthedocs.io/en/latest/) suite (Quinlan and Hall, 2010) was used in read quantification (Figure 1, sub-pipeline 2). The command bedtools coverage generated raw count for TEs while multi-mapping reads matching to n-places (e.g. 10) was recorded n-times (i.e. 10). It also counted number of bases of a TE locus covered by reads (covered bases of TE). Furthermore, this sub-pipeline incorporated bedtools intersect to calculate the number of bases of a read overlapping with an individual TE (mapped bases of read), which were summed for each TE locus to estimate the average read depth of an individual TE’s mapped region (i.e. only the region covered by reads, not the entire annotated feature). This was calculated as follows:
For n reads mapping to a TE locus, and i as integer from 0 to n, f(i) = mapped bases of read i.
Average read depth of an individual TE’s mapped region = <img src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Csum_%7Bi%3D1%7D%5En%20f(i)%20%7D%5Ctext%7BBases%20of%20a%20TE%20locus%20covered%20by%20reads%7D%0A"> .
In order to exclude TE loci that were covered by reads in a sparse and scattered way, a cut-off threshold of the average read depth 5 was adopted in addition to the 10 read count threshold. Examples of the filtering step of this sub-pipeline were illustrated in Figure 2.

<img src="https://github.com/ting-hsuan-chen/TE_ExpressionCandidate/blob/main/Filter_for_bedtools_subpipeline.jpg" width="600">

__Figure 2.__ Examples for the filter step of sub-pipeline 2
<br>
<br>

### Sub-pipeline 3
(Corresponding to scripts 01.02_bwa.bash, 01.05_TEFingerprint.bash, and 01.08_filter_subpipeline_3.r)

The third part (Figure 1, sub-pipeline 3) specifically collects TEs had transcription across the boundaries of the element. The software [TEFingerprint](https://github.com/PlantandFoodResearch/TEFingerprint) (Plant and Food Research) was originally designed for identifying unannotated insertions in genomes using paired-end short fragment DNA sequence data. Here it was applied to capture TE loci internally mapped by multi-mapping reads only yet the read mates, as known as danglers, were uniquely aligned to a location near the insertion site (Figure 3). Although TEFingerprint also has the option to use reads sit across the junction of a TE locus, this function was disabled in sub-pipeline 3 as the utilization of htseq-count in sub-pipeline 1 has covered this scenario. Sub-pipeline 3 applied the standard TEFingerprint pipeline where reads were mapped against the collection of total annotated TE sequences of the investigated species using [BWA](https://github.com/lh3/bwa) (Li and Durbin, 2009). Subsequently the mates of TE-mapped reads were aligned to the reference genome before calculating read count of dangler clusters. Only clusters containing more than 10 dangler reads were kept to test for the intersection of dangler clusters and annotated TEs using bedtools intersect. The candidates need to show more than 10 dangler reads and more than 10 reads mapping internally (counted by bedtools coverage).

<img src="https://github.com/ting-hsuan-chen/TE_ExpressionCandidate/blob/main/TEFingerprint_dangler_reads.jpg" width="600">

__Figure 3.__ TE-related dangler reads counted by TEFingerprint
<br>
<br>

### Summarizing TE expression candidates
(Corresponding to the script 01.09_merge.r)
After excluding TEs did not show enough signal of transcription, potentially expressed TEs from the three sub-pipelines were collected together as a pool of expression candidates.

### Post-pipeline analysis
After identifying expressed TE loci, users can then analyse their characteristics (e.g. integrity or sorting by TE family) or location in relation with annotated genes. Script examples of these analysis can be found in the folder [docx](https://<!>github.com/ting-hsuan-chen/TE_ExpressionCandidate/tree/main/docx).


### References
__Anders S, Pyl PT, Huber W.__ HTSeq--a Python framework to work with high-throughput sequencing data. Bioinformatics. 2015;31:166–9.

__Kim D, Langmead B, Salzberg SL.__ HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12:357–60.

__Li H, Durbin R.__ Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009;25:1754–60.

__Plant and Food Research.__ TEFingerprint, https://github.com/PlantandFoodResearch/TEFingerprint. 2019. 

__Quinlan AR, Hall IM.__ BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26:841–2.

