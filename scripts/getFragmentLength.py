#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:33:08 2023

@author: Marie Van Hecke
"""

import argparse
import HTSeq
from os import listdir
from os import makedirs
from os import stat
from os.path import exists
from multiprocessing import Pool
from multiprocessing import Manager
import tqdm
from functools import partial
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import random

def writeFragLen(outFile, queue):

    out = open (outFile, "w")

    while 1:
        m = queue.get()
        
        if (m == "stop"):
            out.close()
            return

        for fragLen in m:
            out.write(str(fragLen))
            out.write("\n")

def getReference(gtfFile):

    #get exons from gtf file    
    gtf = HTSeq.GFF_Reader(gtfFile)
    nonSplicedGenes = {}
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    
    gene = ""
    numTranscripts = 10
    geneExons = []
    
    for feature in gtf:
        if feature.type == "gene":
            
            #finish off previous gene
            #only save non-spliced genes (= only 1 transcript)
            if numTranscripts == 1:
                nonSplicedGenes[gene] = geneExons
                for iv in geneExons:
                    iv.strand = "."
                    exons[iv] += gene
            
            #reset values for new gene
            gene = feature.name
            numTranscripts = 0
            geneExons = []
        
        if feature.type == "transcript":
            numTranscripts += 1
        
        if feature.type == "exon":
            geneExons.append(feature.iv)
    
    #finish off last gene 
    if numTranscripts == 1:
        nonSplicedGenes[gene] = geneExons
        for iv in geneExons:
            iv.strand = "."
            exons[iv] += gene
                    
    return nonSplicedGenes, exons

def assignReadsToGene(read1, read2, exons):

    iset = None
    minPos = np.inf
    maxPos = 0

    #assign read 1
    for cigop in read1.cigar:
        #skip parts of alignment that are not alignment match, sequence match or sequence mismatch
        if not (cigop.type == "M" or cigop.type == "=" or cigop.type == "X") :
            continue
        minPos = min(minPos, cigop.ref_iv.start)
        maxPos = max(maxPos, cigop.ref_iv.end)
        for iv,step_set in exons[cigop.ref_iv].steps():
            if (len(step_set) > 0):
                if iset is None: #add to empty set
                    iset = step_set.copy()
                else:
                    iset.intersection_update(step_set) #only keep intersection of elements

    #assign read2
    for cigop in read2.cigar:
        #skip parts of alignment that are not alignment match, sequence match or sequence mismatch
        if not (cigop.type == "M" or cigop.type == "=" or cigop.type == "X") :
            continue
        minPos = min(minPos, cigop.ref_iv.start)
        maxPos = max(maxPos, cigop.ref_iv.end)
        for iv,step_set in exons[cigop.ref_iv].steps():
            if (len(step_set) > 0):
                if iset is None: #add to empty set
                    iset = step_set.copy()
                else:
                    iset.intersection_update(step_set) #only keep intersection of elements

    #read does not map to any gene
    if iset==None:
        return {}, minPos, maxPos
    return iset, minPos, maxPos

def readBam(sample, exons, queue):

    #read in bam file to get aligned reads
    try:
        bam_reader = HTSeq.BAM_Reader(sample)
        paired_bam_reader = HTSeq.pair_SAM_alignments_with_buffer(bam_reader)
    except:
        print("file " + sample + " could not be read, check if it exists and if it contains reads.")
        return
    
    ctr = 0
    fragLens = []
    for read1, read2 in paired_bam_reader:
        
        #skip unpaired reads & non primary alignments
        if (read1 is None or read2 is None or read1.not_primary_alignment or read2.not_primary_alignment):
            continue

        try:
            #get set of genes to which read maps
            mappingGenes, minPos, maxPos = assignReadsToGene(read1, read2, exons)
            
            
            if len(mappingGenes) == 1:
                
                #get length of fragment
                gene = list(mappingGenes)[0]
                
                completeInterval = HTSeq.GenomicInterval(read1.iv.chrom, minPos, maxPos, ".")
                fragLen = 0
                for iv, val in exons[completeInterval].steps():
                    if gene in val:
                        fragLen += iv.length
                fragLens.append(fragLen)
                
                pass
            
                ctr += 1
            
        except KeyError:
            #chrom of read not in reference
            pass
    
    queue.put(fragLens)

def removeOverlapping(nonSplicedGenes, exons):
    
    overlappingGenes = set()
    nonOverlappingExons = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    
    #get all overlapping genes
    for iv, val in exons.steps():
        if (len(val) > 1):
            overlappingGenes.update(val)
            
            #remove gene from nonSplicedGenes
            for gene in val:
                if gene in nonSplicedGenes:
                    nonSplicedGenes.pop(gene)
    
    for gene in nonSplicedGenes:
        for iv in nonSplicedGenes[gene]:
            nonOverlappingExons[iv] += gene
    
    return nonSplicedGenes, exons

def _parse_arguments():
    pa = argparse.ArgumentParser(
        add_help=False,
    )
    args, argv = pa.parse_known_args()

    pa = argparse.ArgumentParser(
        parents=[pa],
    )

    pa.add_argument(
        "--bam_dir",
        dest="bamDir",
        type=str,
        help="path/to/bamDir containing the STAR output bam files",
    )

    pa.add_argument(
        "-n",
        "--CPU",
        dest="nCPU",
        type=int,
        default=4,
        help="number of parallel CPU processes (default = 4)",
    )

    pa.add_argument(
        "--gtfFile",
        dest="gtfFile",
        type=str,
        help="path/to/gtfFile",
    )
    
    pa.add_argument(
        "--outDir", 
        dest="outDir",
        type=str,
        help="path/to/outDir",
    )
    
    pa.add_argument(
        "--numBams",
        dest="numBams",
        type=int,
        default=0,
        help="number of bam files to use (default = all files)",
    )
    
    pa.add_argument(
    	"--seed",
    	dest="seed",
    	type=int,
    	default=1234,
    	help="seed for reproducability",
    )

    pa.add_argument(
        "--readLen",
        dest="readLen",
        type=int,
        help="length of 1 read",
    )

    args = pa.parse_args()

    return args

def main(): 
    
    args = _parse_arguments()

    makedirs(args.outDir, exist_ok=True)
    random.seed(args.seed)

    #get non-spliced genes
    print("getting non-spliced genes")
    nonSplicedGenes, exons = getReference(args.gtfFile)
    
    #filter out overlapping genes
    print("filtering overlapping genes")
    filteredGenes, exons = removeOverlapping(nonSplicedGenes, exons)
    
    #get all samples
    bamfiles = [(args.bamDir + "/" + file) for file in listdir(args.bamDir) if file.endswith(".bam")]
    
    #limit number of bamfiles to analyse
    numBams = args.numBams
    if numBams == 0 or numBams > len(bamfiles):
        numBams = len(bamfiles)
    bamfiles = random.sample(bamfiles, numBams)
    
    print("getting fragment lengths for ", numBams, " bams")
    
    manager = Manager()
    writeQueue = manager.Queue()
    
    
    nCPU = min(args.nCPU, len(bamfiles))    #avoid idle threads
    pool = Pool(nCPU + 1)                   #1 additional process for writer
    pool.apply_async(writeFragLen, (args.outDir + "/fragLen" + str(numBams) + ".txt", writeQueue, ))

    #read bams to get fragment length
    for _ in tqdm.tqdm(pool.imap_unordered(partial(readBam, exons = exons, queue = writeQueue), bamfiles), total=len(bamfiles)):
        pass
    
    writeQueue.put("stop")
    
    pool.close()
    pool.join()
    
    print("fitting negative binomial model")
    
    #fit model to fragmentLength, and get mean fragmentlength
    fragLen = pd.read_table(args.outDir + "/fragLen" + str(numBams) + ".txt", header = None)
    fragLen = fragLen[0]

    #select upper and lower limit for values
    upLim = np.quantile(fragLen, 0.99)      #remove top 1% of fragments
    lowLim = args.readLen                   #remove fragments with length < readLen
    fragLen = fragLen[(fragLen < upLim) & (fragLen >= lowLim)]

    #fit negative binomial model
    bounds = [(0, upLim), (0, 1)]
    NBfit = stats.fit(stats.nbinom, fragLen, bounds)
    
    #plot model fit
    x = np.arange(0, upLim)
    y = stats.nbinom.pmf(x, NBfit.params.n, NBfit.params.p)
    plt.hist(fragLen, bins = 200, density = True)
    plt.plot(x, y, label = "modelFit")
    plt.legend()
    plt.xlabel("fragment length (bp)")
    plt.ylabel("density")
    plt.title("fragment length distribution")
    plt.savefig(args.outDir + "/fittedModel" + str(numBams) + ".png")
    
    #get mean and overdispersion
    n = NBfit.params.n
    p = NBfit.params.p
    NBmean = (n*(1-p))/p
    NBdisp = 1/(1-p)
    NBmax = x[y == max(y)][0]
    
    #write params to file
    f = open(args.outDir + "/modelParams" + str(numBams) + ".tsv", "w")
    f.write("n\t" + str(n) + "\n")
    f.write("p\t" + str(p) + "\n")
    f.write("mean\t" + str(NBmean) + "\n")
    f.write("overdispersion\t" + str(NBdisp) + "\n")
    f.write("maxProb\t" + str(NBmax) + "\n")
    f.close()

    #write frequency to file
    fragFreq = pd.Series(fragLen).value_counts().rename_axis('val').reset_index(name='counts')
    fragFreq.to_csv(args.outDir + "/fragLenFreq" + str(numBams) + ".tsv", sep = "\t", index = False)

    
    
    
if __name__ == "__main__":
    main()
