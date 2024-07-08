# ELLIPSIS

A tool for splice graph construction and alternative splicing detection.

## Prerequisites
In order to build, you need the following software:
* conda
* the Eigen library (>3.4) : see https://eigen.tuxfamily.org/
* the BOOST library : see https://www.boost.org/

## Compilation
Clone the GitHub repository:

```bash
git clone https://github.ugent.be/mhecke/ELLIPSIS.git
```

First, create the conda environment using
```bash
conda env create -f ELLIPSIS_environment.yml
conda activate ELLIPSIS_env
```
Always build and run the tool inside this conda environment.

Next, compile the C++ code as follows:
```bash
cd ELLIPSIS
mkdir build
cd build
cmake .. -DEIGEN3_INCLUDE_DIR=/path/to/eigen -DBOOST_INCLUDE_DIR=/path/to/boost
make
```

After compilation, the binary can be found in:
```bash
ELLIPSIS/build/src/ELLIPSIS
```

## Input 

* reference gtf file
* output files from STAR alignment
* nearest neighbors per cell
* raw gene count file

In [scripts](scripts/README.md), an example pipeline is provided as a guideline to obtain these required inputs.
For examples of each of these files, see [testData](testData).

## Usage

Note1: The tool requires reading and writing lots of files simultaneously.
Some systems limit the number of open files.
On linux systems, you can set this limit by using:
    ```bash
    ulimit -n <new limit>
    ```

Note2: It is advised to use STAR in 2pass mode. To avoid removal of novel splice sites detected during the first pass of STAR, use `-novelJunctionFile`.

Note3: If you are using single end reads, make sure to use the `-singleEnd` option.


Running ELLIPSIS to get robust PSI values:

```bash
ELLIPSIS -gtf </path/to/ref.gtf> \
         -alignmentDir </path/to/alignment>  \
         -neighborFile </path/to/neighbors.tsv> \
         -countFile </path/to/scCounts.tsv> \
         -readLength <RL> \
         [options] 
```

required arguments: 
* `/path/to/ref.gtf` the gtf file (from Ensembl) containing the reference genome annotation <br />
* `/path/to/alignment` the directory containing the STAR output files (bam files and SJ.out.tab files) <br />
* `/path/to/neighbors.tsv` the file containing a list of k nearest neighbors (tab separated)
* `/path/to/countfile` the file containing the geneCounts (tab separated) <br />
* `RL` the average length of a read (for PE reads: length of 1 mate) <br />

 \[options without argument\]: 

* `-verbose`: create additional log files during PSI calculations (logPSI and logAlpha) <br/>
* `-singleEnd`: single end reads [default PE reads] <br />

 \[options with 1 argument\]: 

* `-novelJunctionFile` : file containing novel junctions added to the STAR reference prior to 2nd pass of STAR <br />
* `-minJunctionCount` : minimum number of junction spanning reads for 1 cell to be considered expressed [default = 5] <br />
* `-minJunctionCells` : minimum number of cells for which non annotated junction is expressed for it to be accepted [default = 10] <br />
* `-minDepth` : minimum read depth required to consider a cell for computing PSI values [default = 10]  <br />
* `-wObs` : weight assigned to observed depths [default = 1] <br />
* `-wSrcSink` : weight assigned to impose src/sink have PSI = 100% [default = 1]  <br />
* `-wFlow` : weight assigned to impose conservation of flow of PSI values [default = 6]  <br />
* `-wClust` : weight assigned to favor intra-cell type similarity [default = 2]  <br />
* `-maxIter` : maximum number of iterations in compute PSI step [default = 100] <br />
* `-maxLowQual` : maximum fraction of low quality read mappings to gene [default = 0.1] <br />
* `-maxPaths` : maximum number of source to sink paths to compute PSI values (use inf to remove maxPaths filter) [default = 1e9] <br />
* `-chunkSize` : maximum number of files read at once for SJMerge [default = 1000] <br />
* `-outDir` : output directory [default = this directory] <br />
* `-run` : run name, use when you are performing several runs<br />
* `-CPU` : number of CPU threads [default = 4] <br />
* `-selectGenes`: file containing genes to consider. If not provided, PSI values are computed for all genes possible <br />


* `-trueClusterFlow` : directory containing trueClusterFlow files for each gene <br/>
* `-clusterFile` : clusters.tsv the file containing a list of cells and their corresponding cluster (tab separated) <br/>
* `-trueGTF`: gtf file containing only transcripts that are present, if option not given, no trueGraphs are build <br/>

The last 3 options can be used for simulated data, where the ground truth is known.

The program consists of a few steps:
When rerunning the program, the existence of the output directories is checked. 
If they exist, the according step is skipped. 
In order to force re-running these steps, remove the according directories from the output directory.
Of you are performing several runs with slightly different parameters that only affect the computePSI step, you can use the `-run` parameter. 
This will put all the outputs of the computePSI in a separate directory, without needing to remove the previous outputs of that step.

### 1. BuildRefGraphs
Builds splice graphs for each gene from the reference genome (gtf file).
The results are stored in the `annotGraphs` directory.

For each gene, a file is added with the first line containing gene information,
followed by each junction present in the graph. <br />
<pre>
chromosome   geneID   geneName   strand[+/-]   geneStart   geneEnd
E1start-E1end -> E2start-E2end
E2start-E2end -> E3start-E3end
...
</pre>

Exon positions are always ordered using Estart < Eend, even for genes on the reverse strand.
All positions are 1-based.

### 2. MergeSJFiles
STAR outputs observed annotated and novel splice junctions per cell. 
In this step, the reported junctions are filtered to only obtain novel splice junctions.
Additionally, the reported novel junctions are aggregated per gene, and the number of uniquely mapping reads are stored per cell. 
The resulting count matrices are written to the `novelSJCounts` directory.
Genes without any observed novel junctions are not stored.

Additionally, we also report the junctions that could not be attributed to a gene, and are therefore ignored.
They can be found in junctionsWithoutGene.tsv. 
The junctions that could not be attributed to 1 gene, but were added to multiple genes are listed in junctionsMultipleGenes.tsv.
These junctions are added to the novel junctions of each gene.

The genomic positions of the junctions use the same format as reported in the STAR SJ.out.tab files.
The junctions are named as follows: chromosome_startPos_endPos, with:
* startPos : first base of the skipped intron (1-based)
* endPos : last base of the skipped intron (1-based)

The startPos and endPos are always ordered using startPos < endPos, even for genes on the reverse strand.

### 3. AddDeNovo
Adds de novo splice events to the graphs.
If an observed un-annotated junction occurs more than `minJunctionCount` in more than `minJunctionCells` samples, the corresponding junction is added to the splice graph of the corresponding gene.
Novel junctions that are entirely contained in any gene, are omitted.
Novel junctions that map are contained in multiple genes are added to each of the corresponding splice graphs.

The results are stored in the `deNovoGraphs` directory, using the same file format as described in [BuildRefGraphs](#1-buildrefgraphs).

### 4. GetDepth
In this step, gene level counts are computed, as well as exon-specific and junction-specific depths.
Reads that map ambiguously to multiple genes, will be proportionally assigned to each gene given the gene expression.

The junction counts are stored per gene in the `allJunctionCounts` directory.
Using the following format: 
<pre>
        J1Low_J1High   J2Low_J2High     J3Low_J3High
cell1.bam   5              0                7
cell2.bam   3              8                0
</pre>
Each number in the count table represents a read in the cell mapping a junction.

The depth per exon is stored per gene in the `exonDepth` directory.
Using the following format:
<pre>
        E1Low-E1High   E2Low-E2High     E3Low-E3High
cell1.bam   5              0                7
cell2.bam   3              8                0
</pre>
The numbers in this count table represent the sum of the per base exon coverage.
Therefore, larger exons will typically have higher exonDepths in this table.

We also get gene level counts stored in the `geneCounts` directory.
Using the following format:
<pre>
cell    counts
cell1.bam   10           
cell2.bam   300           
</pre>
Where each number corresponds to the number of reads mapping to that gene.

Finally, the fraction of low quality read alignments that map to each gene is stored in the `lowQual` directory.
<pre>
cell    lowQualFrac
cell1.bam   0.01
cell2.bam   0.15
</pre>


### 5. SimplifyGraphs

Removing unobserved exons and junctions from the graph, while making sure the simplified graph remains connected. 
Exons and junctions that did not occur in the `exonDepth` resp. `allJunctionCounts` file, are removed.

In case the `-selectGenes` option was given, we only consider the genes in the given genelist.

The resulting simplified graphs are stored in the `simplifiedGraphs` directory using the same file format as described in [BuildRefGraphs](#1-buildrefgraphs).

### 6. computePSI
The PSI values are estimated based on the local read coverage, conservation of flow of PSI values in the graphs, and intra-cell type similarity.

The resulting PSI values are reported in the `PSI` directory, using the following format.
<pre>
cell        E1      E2      ...     J1      J2
cell1.bam   1       0.52            0.53    0.2
cell2.bam   1       0.49            0.51    0.23
</pre>

### 7. \[optional\] Get ground truth
This step is only performed if trueGTF and trueClusterFlow are given.
Here, the most simple splicegraphs representing all splice variation in the ground truth are build, and stored in the `trueGraphs` directory.
The splicegraph files follow the same structure as described in [BuildRefGraphs](#1-buildrefgraphs)

Additionally, also the true PSI values for each individual cell are computed, and stored in the `trueCellPSI`directory.

### 8. \[optional\] Get error
This step is only performed if trueGTF and trueClusterFlow are given.
The average PSI error and the average absolute PSI error made per gene are computed.

