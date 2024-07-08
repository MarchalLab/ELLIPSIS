#!/bin/bash

N=1 #default value for number of cores
starIdx=None
while getopts i:s:f:g:o:n: flag
do
    case "${flag}" in
        i) inDir=${OPTARG};;
        s) starIdx=${OPTARG};;
        f) fastaFile=${OPTARG};;
        g) gtfFile=${OPTARG};;
        o) outDir=${OPTARG};;
        n) N=${OPTARG};;
    esac
done


if [ ! -d "$inDir" ]; then
    echo "$inDir does not exist."
    exit
fi

## create STAR index if it does not yet exist
if [ $starIdx == None ] | [ ! -d $starIdx ]; then
    
    echo "generating STAR index"
    
    if [ ! -f "$fastaFile" ]; then 
        echo "$fastaFile does not exist."
        exit
    fi
    
    if [ ! -f "$gtfFile" ]; then 
        echo "$gtfFile does not exist."
        exit
    fi
    
    #create STAR index for reference
    STAR --runThreadN $N \
        --runMode genomeGenerate \
        --genomeDir $outDir/STAR_index \
        --genomeFastaFiles $fastaFile \
        --sjdbGTFfile $gtfFile
    
    starIdx=$outDir/STAR_index
fi

## 1st pass
mkdir $outDir/pass1
STAR --genomeLoad LoadAndExit --genomeDir $starIdx
for pattern in _1.fa _1.fa.gz _1.fasta _1.fasta.gz _R1.fa _R1.fa.gz _R1.fasta _R1.fasta.gz _1.fq _1.fq.gz _1.fastq _1.fastq.gz _R1.fq _R1.fq.gz _R1.fastq _R1.fastq.gz
do
    pattern2=${pattern//1/2}
    
    for file in ${inDir}/*${pattern}
    do
        if [ -f $file ]
        then
            #get samplename
            sample=${file##*/}
            sample=${sample%${pattern}}

            if [[ $pattern == *.gz ]]
            then 
                STAR --runThreadN $N \
                     --genomeLoad LoadAndKeep \
                     --genomeDir $starIdx \
                     --readFilesIn ${inDir}/${sample}${pattern} ${inDir}/${sample}${pattern2} \
                     --outFileNamePrefix $outDir/pass1/${sample} \
                     --readFilesCommand zcat \
                     --outSAMtype None
            else
                STAR --runThreadN $N \
                     --genomeLoad LoadAndKeep \
                     --genomeDir $starIdx \
                     --readFilesIn ${inDir}/${sample}${pattern} ${inDir}/${sample}${pattern2} \
                     --outFileNamePrefix $outDir/pass1/${sample} \
                     --outSAMtype None
            fi
        fi
    done
    
done
STAR --genomeLoad Remove --genomeDir $starIdx

#filtering junctions
echo "Number of new junctions before filtering"
cat $outDir/pass1/*SJ.out.tab | awk '($6==0)' | cut -f1-6 | sort | uniq |wc -l

echo "Number of new junctions with a non-canonical intron motif"
cat $outDir/pass1/*SJ.out.tab | awk '($6==0 && $5 > 0)' | cut -f1-6 | sort | uniq |wc -l

echo "Number of new junctions without enough read support (<5 in all cells)"
cat $outDir/pass1/*SJ.out.tab | awk '($6==0 && $7>5)' | cut -f1-6 | sort | uniq | wc -l

#remove junctions with
 #non cannonical intron motif ($5>0)
 #junctions that are already annotated ($6==0) -> these are already in gtf file
 #not enough uniquely mapping reads ($7>5) -> limit holds for each file separately, but if enough reads in 1 cell, the junction will be added for all cells in the dataset
cat $outDir/pass1/*SJ.out.tab | awk '($5 > 0 && $6==0 && $7>5)' | cut -f1-6 | sort | uniq > $outDir/pass1/SJ.filtered.tab

echo "Number of new junctions after filtering"
wc -l < $outDir/pass1/SJ.filtered.tab

#create updated STAR reference
STAR --runThreadN $N \
     --runMode genomeGenerate \
     --genomeDir ${outDir}/STAR_index2 \
     --genomeFastaFiles $fastaFile \
     --sjdbGTFfile $gtfFile \
     --sjdbFileChrStartEnd $outDir/pass1/SJ.filtered.tab

#2nd pass
mkdir $outDir/alignment
STAR --genomeLoad LoadAndExit --genomeDir ${outDir}/STAR_index2
for pattern in _1.fa _1.fa.gz _1.fasta _1.fasta.gz _R1.fa _R1.fa.gz _R1.fasta _R1.fasta.gz _1.fq _1.fq.gz _1.fastq _1.fastq.gz _R1.fq _R1.fq.gz _R1.fastq _R1.fastq.gz
do
    pattern2=${pattern//1/2}
    
    for file in ${inDir}/*${pattern}
    do
        if [ -f $file ]
        then
            #get samplename
            sample=${file##*/}
            sample=${sample%${pattern}}

            if [[ $pattern == *.gz ]]
            then 
                STAR --runThreadN $N \
                     --genomeLoad LoadAndKeep \
                     --genomeDir $outDir/STAR_index2 \
                     --readFilesIn ${inDir}/${sample}${pattern} ${inDir}/${sample}${pattern2} \
                     --outFileNamePrefix ${outDir}/alignment/${sample} \
                     --readFilesCommand zcat \
                     --outFilterType BySJout \
                     --limitBAMsortRAM 10000000000 \
                     --outSAMtype BAM SortedByCoordinate
            else
                STAR --runThreadN $N \
                     --genomeLoad LoadAndKeep \
                     --genomeDir $outDir/STAR_index2 \
                     --readFilesIn ${inDir}/${sample}${pattern} ${inDir}/${sample}${pattern2} \
                     --outFileNamePrefix ${outDir}/alignment/${sample} \
                     --outFilterType BySJout \
                     --limitBAMsortRAM 10000000000 \
                     --outSAMtype BAM SortedByCoordinate
            fi
        fi
    done
    
done
STAR --genomeLoad Remove --genomeDir ${outDir}/STAR_index2

#index file
for file in $outDir/alignment/*Aligned.sortedByCoord.out.bam 
do

    #get samplename
    sample=${file##*/}
    sample=${sample%Aligned.sortedByCoord.out.bam}
    
    (
        #index bam file
        samtools index $file
        
    )&
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then 
    wait -n; fi
done
wait

## filtering (remove non-primary alignments & unmapped reads)
mkdir $outDir/alignment_filtered
for file in $outDir/alignment/*Aligned.sortedByCoord.out.bam 
do
        (
        
        #get samplename
        sample=${file##*/}
        sample=${sample%Aligned.sortedByCoord.out.bam}
    
        #filter out non-primary alignments (256) and unmapped reads (4), and only include reads with proper pair (2)
        samtools view -h -b -F 256 -F 4 -f 2 $file > $outDir/alignment_filtered/${sample}.bam

        #index bamfile
        samtools index $outDir/alignment_filtered/${sample}.bam 
        
        )&
        
        if [[ $(jobs -r -p | wc -l) -ge $N ]]
        then 
            wait -n
        fi

done
wait
