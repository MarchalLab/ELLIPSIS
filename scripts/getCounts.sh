#!/bin/bash

while getopts i:o:g:n: flag
do
    case "${flag}" in
        i) bamDir=${OPTARG};;
        o) outPath=${OPTARG};;
        g) gtf=${OPTARG};;
        n) N=${OPTARG};;
    esac
done

if [ ! -d "$bamDir" ]; then
    echo "$bamDir does not exist."
    exit
fi

if [ ! -f "$gtf" ]; then
    echo "gtf file $gtf does not exist."
    exit
fi

thisDir=$(dirname "$0")

## counts
mkdir $outPath

for file in $bamDir/*.bam 
do

    #get samplename
    sample=${file##*/}
    
    (
        #gene level counts using HTSeq
        htseq-count \
                -m intersection-nonempty \
                -s no \
                -r pos \
                $file \
                $gtf \
                -q \
                > $outPath/${sample}.counts
        
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then 
    wait -n; fi
done
wait

echo "merging counts"

#merge HTSeq counts into 1 file
Rscript $thisDir/mergeCounts.R -d $outPath


