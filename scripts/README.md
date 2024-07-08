# preprocessing data

Starting from the raw fastq files of a single cell dataset, we present a pipeline that can be used in order to generate all the required input data for ELLIPSIS.
The full pipeline can be found in `/ELLIPSIS/scripts/pipeline/preprocessing.sh`.


## prerequisites
For most steps in the preprocessing pipeline, dedicated conda environments are required. 
These can be created using the following commands:
```bash
conda env create -f bedtools_environment.yml
conda env create -f scDataProcessing_environment.yml
```

## required inputs
* `inDir` containing paired end read files
* `gtfFile` the reference gtf file (e.g. from ensembl)
* `fastaFile` the reference fasta file (e.g. from ensembl)
* `outDir` the path where all output is written

## filter readthrough transcripts from gtf file
It is advised to use a gtf file that does not contain readthrough transcripts. 
Therefore, we filter out the readthrough transcripts from the gtf file.

## STAR alignment
The raw reads are aligned using STAR in 2pass mode. 
For computational reasons, we filter out non-primary alignments, unmapped reads, and reads without proper pair, as these are not considered for further computations.
Reads are sorted by coordinate to be compatible with HTSeq-count.

## gene level counts
We also require a single count file containing gene level counts for each cell.
In this pipeline, we use HTSeq-count to obtain gene level counts.
Other RNA count tools can be used, but the resulting file should be formatted as follows:
  <pre>
          cell1.bam   cell2.bam     cell3.bam
  gene1   5              0                7
  gene2   3              8                0
  </pre>

## Nearest neighbors
We require a list of neighboring cells for each cell in the dataset.
The neighboring cells can be computed based on the raw gene counts using Seurat.
Other approaches can be used to find neighboring cells: e.g. selecting all cells that are from the same cell type as each others neighbors.
The resulting neighbor file should be formatted as follows (tab-separated):
<pre>
cell1.bam   neighbor1.bam   neighbor2.bam   neighbor3.bam
cell2.bam   neighbor1.bam   neighbor2.bam   neighbor3.bam
</pre>

For examples of the input file, see [testData](/testData)






