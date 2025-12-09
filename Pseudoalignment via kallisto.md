# Pseudoalignment via kallisto

## Table of Contents

- Brief Overview of RNA Sequencing
- Downside of Traditional Alignment
- Pseudoalignment
- Traditional Alignment Tools vs. Kallisto
- Pseudoalignment Costs and Benefits
- Kallisto Tutorial
- Conclusion
- References
## Brief Overview of RNA Sequencing

&nbsp;&nbsp;&nbsp;&nbsp;RNA sequencing data can reveal important information about what genes are expressed in organisms at specific locations, times, and quantities. In order to decode the transcriptome, RNA molecules must be dissociated from the samples from which they were taken. Then, complementary DNA is sequenced and fragmented from the RNA molecules to produce a sequencing library.  Once the DNA is sequenced, the resultant reads can be mapped to a reference genome or transcriptome to determine specific gene expression. RNA sequencing is particularly useful for characterizing cell and tissue types, disease states, and more (Jiang et al).

## Downsides of Traditional Alignment

&nbsp;&nbsp;&nbsp;&nbsp; The RNA-sequencing data analysis pipeline involves the alignment of sequenced reads to a reference genome or transcriptome, the quantification of the number of reads aligned to each gene in the reference, the identification of differentially expressed genes, and the discernment of enriched biological pathways. 

&nbsp;&nbsp;&nbsp;&nbsp; The alignment and quantification steps of RNA-sequencing data analysis are incredibly time and resource intensive — traditional alignment tools such as TopHat and Cufflinks can take tens of hours and require significant computing power in order to complete alignment and quantification. This is because traditional alignment requires the matching of each individual read to a particular location on a reference genome, which is computationally expensive when millions of reads are involved. As datasets become progressively larger and analysis takes longer, it becomes imperative to reduce processing time for quantification.

&nbsp;&nbsp;&nbsp;&nbsp; Some alternative methods of quantification have been introduced, but they increase time efficiency at the expense of accuracy. One such method called eXpress involves a streaming algorithm that calculates the probabilities for multiple sites that a particular read can map to based on prior approximations. Another tool, htseq-count, only retains a subset of reads that clearly map to only one site on a reference. Additionally, the tool Sailfish utilizes hashing and lookup tables to match k-mers from a dataset, which obscures essential information that can only be discerned from full-length reads.

&nbsp;&nbsp;&nbsp;&nbsp; In response to the challenges associated with alignment and quantification, the tool “kallisto” was created to maximize time efficiency without sacrificing accuracy by employing a technique called pseudoalignment (Bray et al).

## Pseudoalignment

&nbsp;&nbsp;&nbsp;&nbsp; Unlike traditional alignment tools, Kallisto uses a “pseudoalignment” approach that focuses on identifying which transcripts a read is compatible with rather than determining the precise transcript alignment position. By bypassing the need for full alignment of reads to a genome or transcriptome, pseudoalignment greatly reduces processing time. 

&nbsp;&nbsp;&nbsp;&nbsp; Kallisto starts by constructing an index of the reference transcriptome, using a FASTA file as input and decomposing all transcript sequences into smaller segments known as k-mers. These k-mers are connected based on overlaps of k-1 nucleotides, forming a De Bruin graph (T-DBG) in which each k-mer represents a node and overlaps define the edges between them. Linear stretches in the graph correspond to continuous regions of transcripts known as contigs. Once the T-DBG is constructed, kallisto creates its index by storing a hash table that maps each k-mer to the contig it belongs to, along with its position within the contig. 

&nbsp;&nbsp;&nbsp;&nbsp; To begin pseudoalignment for sequenced samples, kallisto first decomposes each read into overlapping k-mers. Each k-mer is then looked up in the hash table created from the T-DBG in order to identify its location in the graph, forming a path that traces the read across the T-DBG. From this matching, k-compatibility classes are generated for each k-mer, representing the sets of transcripts that each k-mer could have originated from. For each read, taking the intersection of all its k-mer’s k-compatibility classes then determines the read’s equivalence class, producing the set of transcripts compatible with the entire read. Many reads may share the same equivalence class, allowing kallisto to group up these reads together and count them collectively rather than tracking each individual read. 

&nbsp;&nbsp;&nbsp;&nbsp; After forming these equivalence classes, kallisto uses the Expectation-Maximization (EM) algorithm to maximize the probabilities of selecting reads from specific transcripts, essentially determining the most probable distribution of reads across transcripts (Bray et al).

## Figures
[include comparison, pseudoalignment pictures here]
## Traditional Alignment Tools vs. Kallisto
&nbsp;&nbsp;&nbsp;&nbsp; Traditional alignment-based tools and pseudoalignment tools each offer distinct strengths depending on the downstream analysis goals and characteristics of the experimental data. STAR, one of the most widely used traditional alignment tools, performs full read alignment and is splice aware. Because STAR can map reads across exon–exon junctions, the tool is well suited for analyses involving alternative splicing, datasets with high sequencing depth, and datasets where genes are expressed at low levels. However, STAR’s full alignment process is computationally intensive and significantly slower than the pseudoalignment approach.

&nbsp;&nbsp;&nbsp;&nbsp; In comparison, Kallisto is extremely fast and computationally efficient while still offering expression estimates that correlate strongly with results from traditional aligners. Kallisto’s pseudoalignment strategy makes it a practical choice for large datasets or limited computational resources. Although it is not splice aware, Kallisto performs well when working with samples that have high isoform diversity, datasets with low sequencing depth, or well-annotated, complete transcriptomes without extensive splicing events. Overall, STAR offers advantages when accuracy across exon junctions, full read alignments, or detailed downstream analyses are essential, while Kallisto excels in speed, efficiency, and high quality transcript quantification are priorities (Srivastava and Jiang et al).

## Pseudoalignment Costs and Benefits
&nbsp;&nbsp;&nbsp;&nbsp; Pseudoalignment may be preferable over traditional alignment when the input data involves longer genes that code for proteins. Pseudoalignment and quantification through kallisto has a similar accuracy rate to traditional alignment tools such as HISAT2, but is significantly more time and resource efficient.

&nbsp;&nbsp;&nbsp;&nbsp; However, there are still certain situations where it would be preferable to use traditional alignment methods for quantification. Kallisto underperforms relative to traditional alignment tools when the input data includes low-abundance RNA and/or smaller-length RNA. This may be due to the presence of RNA bases modified after transcription that cannot be mapped to a reference, which then inhibits k-mer counting and mapping (Wu et al and Yi et al).

## Kallisto Tutorial, adapted from the [Pachter Lab](https://pachterlab.github.io/kallisto/starting) 

### Step 1: environment setup
- ``conda activate beng_183``
- ``mkdir kallisto_tutorial``
- ``cd kallisto_tutorial``
- Install kallisto using conda following this [guide](https://pachterlab.github.io/kallisto/download)
[manual picture here]
### Step 2: Building an index 
- Download GTF (comprehensive gene annotation CHR) and FASTA transcripts (all) files from [GENCODE](https://www.gencodegenes.org/human/)
- ``kallisto index -i gencode.v49.transcripts.idx gencode.v49.transcripts.fa.gz``
 [tdbg picture here]
- Downloaded paired end RNA-seq reads from [Run SRR493366](https://www.ebi.ac.uk/ena/browser/view/SRR493366)
[reads here]

### Step 3: Quantification
- ``kallisto quant -i gencode.v49.transcripts.idx -o output -b 100 SRR493366_1.fastq.gz SRR493366_2.fastq.gz ``
NOTE: bootstrap value of 100 is a default value used in examples to estimate the technical variability and uncertainty in transcript abundance quantification but can take much longer to run!
- The results of this kallisto run are placed in the output directory specified after -o. Here are the contents of our output directory after quantifying abundance:
[output here]

### Step 3: Quantification Visualization (Optional)
- To visualize pseudoalignments kallisto needs to be run with the [--genomebam](https://pachterlab.github.io/kallisto/manual) option with a gene annotation file and an optional, but recommended text file containing the length of each chromosome.
- ``kallisto quant -i gencode.v49.transcripts.idx -o kallisto_out -b 30 --genomebam gencode.v49.transcripts.gtf.gz --chromosomes chrom.txt SRR493366_1.fastq.gz SRR493366_2.fastq.gz ``
- --genomebam projects the pseudo alignments to genome coordinates and outputs a bam and bam.bai file in addition to output files in the basic quantification step. The projected pseudoalignments can then be viewed using genome browsers like IGV. 
## Step 4: Analyze abundance.tsv with Sleuth For Multiple Samples (Optional)

- ``conda install --channel bioconda r-sleuth``
- launch R
``create an auxiliary table that describes the experimental design and the relationship between the kallisto directories and the samples. Our example:``

[Insert picture here]

For full guide on using sleuth:
[Getting started with sleuth](https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html)

## Conclusion:

Kallisto is valuable for the proper processing of RNA-sequencing data due to its fast and accurate ability to quantify transcript abundance. Kallisto has assisted with the construction of detailed transcriptomic profiles for various different organisms, and is also used in wrapper packages for rapid RNA-seq analysis. 

## Works Cited

Bray, Nicolas L, et al. “Near-Optimal Probabilistic RNA-Seq Quantification.” Nature Biotechnology, vol. 34, no. 5, 4 Apr. 2016, pp. 525–527, https://doi.org/10.1038/nbt.3519.
Download.” Github.io, 4 Oct. 2019, pachterlab.github.io/kallisto/download. Accessed 17 Nov. 2025.
“Getting Started.” Github.io, 2018, pachterlab.github.io/kallisto/starting.
Jiang, Gao, et al. “A Comprehensive Workflow for Optimizing RNA-Seq Data Analysis.” BMC Genomics, vol. 25, no. 1, 24 June 2024, https://doi.org/10.1186/s12864-024-10414-y.
“Kallisto Manual.” Pachterlab.github.io, pachterlab.github.io/kallisto/manual. 
Yi, Lynn, et al. “A direct comparison of genome alignment and transcriptome pseudoalignment.” bioRxiv 444620; doi: https://doi.org/10.1101/44620 
Pachter Lab. “Apps.” Github.io, 2018, pachterlab.github.io/kallisto/apps. Accessed 17 Nov. 2025.
Pimentel, Harold, et al. “Getting Started with Sleuth.” Github.io, 2025, pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html. Accessed 17 Nov. 2025.
Srivastava, Kriti. “Kallisto vs. STAR: Alignment and Quantification of Bulk RNA-Seq Data.” Www.elucidata.io, 29 May 2023, www.elucidata.io/blog/kallisto-vs-star-alignment-and-quantification-of-bulk-rna-seq-data.   Accessed 16 Nov. 2025.
Wu, Douglas C., et al. “Limitations of Alignment-Free Tools in Total RNA-Seq Quantification.” BMC Genomics, vol. 19, no. 1, 3 July 2018, https://doi.org/10.1186/s12864-018-4869-5. 











