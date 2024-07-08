#!/bin/bash

while getopts i:o: flag
do
    case "${flag}" in
        i) inGTF=${OPTARG};;
        o) outPath=${OPTARG};;
    esac
done

if [ ! -f "$inGTF" ]; then
    echo "$inGTF does not exist."
    exit
fi

#get readthrough transcripts
grep readthrough_transcript $inGTF | gtf2bed | cut -f4 > $outPath/readthroughTranscripts.txt

#get unique geneNames of readthrough transcripts
sort $outPath/readthroughTranscripts.txt | uniq > $outPath/readthroughGenes.txt

echo "number of readthrough genes to remove"
wc -l $outPath/readthroughGenes.txt

#remove readthrough transcripts from genome
cp $inGTF $outPath/noReadThrough.gtf
while read gene
do
    grep -v $gene $outPath/noReadThrough.gtf > $outPath/tempFile.gtf
    mv $outPath/tempFile.gtf $outPath/noReadThrough.gtf
done < $outPath/readthroughGenes.txt



