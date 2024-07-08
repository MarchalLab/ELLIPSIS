#!/bin/bash

################################ Input parameters ################################
readDir=path/to/rawReads
gtfFile=path/to/ref.gtf
fastaFile=path/to/ref.fasta
outDir=path/to/out
N=1
##################################################################################

bamDir=$outDir/alignment_filtered
countDir=$outDir/counts
countFile=$outDir/counts/scCounts.tsv

mkdir $outDir
thisDir=$(dirname "$0")


################## filter readthrough transcripts from gtf file ##################
if true
then
    mkdir $outDir/refNoReadThrough
    conda activate bedtools 
    bash $thisDir/removeReadthrough.sh -i $gtfFile -o $outDir/refNoReadThrough
    conda deactivate

    gtfFile=$outDir/refNoReadThrough/noReadThrough.gtf
fi
##################################################################################

############################ align reads using STAR ##############################
if true
then
    conda activate scDataProcessing 
    bash $thisDir/STARalign.sh -i $readDir -f $fastaFile -g $gtfFile -o $outDir -n 4
    conda deactivate
    bamDir=$outDir/alignment_filtered
fi
##################################################################################

############################### gene level counts ################################
if true
then
    conda activate scDataProcessing
    bash $thisDir/getCounts.sh -i $bamDir -o $outDir/counts -g $gtfFile -n $N
    conda deactivate
    countFile=$outDir/counts/scCounts.tsv
fi
##################################################################################

############################### nearest neighbors ################################
if true
then
    Rscript $thisDir/findKNN.R --countFile $countFile --gtfFile $gtfFile --outFile $outFile/neighbors.tsv --numNB 20
fi
##################################################################################
