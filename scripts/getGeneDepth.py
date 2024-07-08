#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =========================================================================================================
# """
# Script to assign aligned reads to genes and get corresponding exon depth, junction depth and gene counts
# """
# =========================================================================================================

import HTSeq
from os import listdir
from os import makedirs
from os.path import exists
from multiprocessing import Pool
from multiprocessing import Manager
import tqdm
import sys
import argparse
import os.path
import traceback
from functools import partial
import pandas as pd

#Helper classes
class CountDict:

    def __init__(self):
        self.counter = {}

    def __iter__(self):
        return iter(self.counter)

    def add(self, key, val):
        if key not in self.counter:
            self.counter[key] = val
        else:
            self.counter[key] += val

    def get(self, key):
        retVal = self.counter.get(key)
        if retVal == None:
            return 0 #if key does not exist: 0 counted
        return retVal

class GeneInfo:

    def __init__(self, line):
        self.chr, self.geneID, self.geneName, self.strand, self.start, self.stop, self.longestTranscriptLen, self.shortestTranscriptLen = line.split('\t')
        self.observableJunctions = set()
        self.consecutiveJunctions = set()
        self.observableExons = set()
        self.cellDepth = {}

    def addExon(self, exon):
        if (exon != "source" and exon != "sink"):
            self.observableExons.add(exon)

    def addJunction(self, line):

        #remove '\n' from end of string
        line = line.strip()

        #get start and stop exon for this edge
        e1, e2 = line.split(" -> ")

        #add both exons to observable exons
        self.addExon(e1)
        self.addExon(e2)

        #ignore non-observable junctions
        if (e1 == "source" or e1 == "sink" or e2 == "source" or e2 == "sink"):
            return

        #get splice junction
        if (self.strand == "+"): #forward strand
            tmp, j1 = e1.split("-")
            j2,tmp = e2.split("-")
        else: #reverse strand
            j2, tmp = e1.split("-")
            tmp, j1 = e2.split("-")

        #convert from 1-based (STAR) to 0-based (HTSeq)
        j1 = int(j1) - 1
        j2 = int(j2) - 1

        #add to junction sets
        if j2-j1 > 1:
            self.observableJunctions.add((j1,j2))
        else:
            self.consecutiveJunctions.add((j1,j2))

    def addCellDepth(self, cell, depth):
        self.cellDepth[cell] = depth

    def getDepth(self, cell):
        return self.cellDepth.get(cell, 0)

    def getOrderedExons(self):
        return sorted(self.observableExons)

    def getOrderedJunctions(self):
        return sorted(self.observableJunctions) + sorted(self.consecutiveJunctions)

class GeneCov:

    def __init__(self, geneInfo):
        self.exonCov = HTSeq.GenomicArray("auto", stranded = False, typecode = "d")
        self.junctionCounts = CountDict()
        self.consecutiveJunctions = geneInfo.consecutiveJunctions
        self.observableJunctions = geneInfo.observableJunctions
        self.highQualCounts = 0
        self.lowQualCounts = 0

    def addSingle(self, read, depth):
        #skip low quality reads
        if (read.aQual < 10):
            self.lowQualCounts += depth
            return

        #get covered positions & junctions
        posCovered = HTSeq.GenomicArray("auto", stranded = False, typecode = "b")
        junctionsCovered = set()
        self.addSingleRead(read, posCovered, junctionsCovered)

        #add observed coverage
        for iv, val in posCovered.steps():
            if (val):
                self.exonCov[iv] += depth

        #add observed junctions
        for junction in junctionsCovered:
            self.junctionCounts.add(junction, depth)

        #add raw gene counts
        self.highQualCounts += depth

    def addPaired(self, read1, read2, depth):

        #skip low quality reads
        if (read1.aQual < 10 or read2.aQual < 10):
            self.lowQualCounts += depth
            return

        #avoid overlapping mate pair to be counted twice
        posCovered = HTSeq.GenomicArray("auto", stranded = False, typecode = "b")
        junctionsCovered = set()
        self.addSingleRead(read1, posCovered, junctionsCovered)
        self.addSingleRead(read2, posCovered, junctionsCovered)

        #add observed coverage
        for iv, val in posCovered.steps():
            if (val):
                self.exonCov[iv] += depth

        #add observed junctions
        for junction in junctionsCovered:
            self.junctionCounts.add(junction, depth)

        #add raw gene counts
        self.highQualCounts += depth

    def addSingleRead(self, read, posCovered, junctionsCovered):

        prev = None
        for cigop in read.cigar:

            #skip parts of alignment that are not alignment match, sequence match or sequence mismatch
            if not (cigop.type == "M" or cigop.type == "=" or cigop.type == "X") :
                continue

            #genome position of (partial) alignment
            cigopiv = cigop.ref_iv
            cigopiv.strand = '.'

            #add count to coverage vector for matches
            posCovered[cigopiv] = True

            #check if read covers artificial junctions
            for p1,p2 in self.consecutiveJunctions:
                if (cigopiv.contains(HTSeq.GenomicPosition(read.iv.chrom, int(p1)))
                       and cigopiv.contains(HTSeq.GenomicPosition(read.iv.chrom, int(p2)))):
                    junctionsCovered.add((p1,p2))

            #check if read covers non-artificial junction
            if prev is not None:
                p1, p2 = prev.end-1, cigopiv.start
                if (p1,p2) in self.observableJunctions:
                    junctionsCovered.add((p1,p2))

            prev = cigop.ref_iv

    def getLowQualFrac(self):
        return self.lowQualCounts / (self.lowQualCounts + self.highQualCounts)

#write to files asynchronous
def writeFiles(outDir, queue, isPauseQueue):

    makedirs(outDir, exist_ok=True)

    while 1:
        m = queue.get()
        if (m == "stop"):
            return

        fileName, header, line = m

        if isPauseQueue:
            #start reading new bams if less than 10 000 messages remaining
            if ((not unpaused.is_set()) and queue.qsize() < 10000):
                unpaused.set()

            #pause reading new bams if queue is too large (more than 1 Mi messages)
            if (unpaused.is_set() and queue.qsize() > 1000000):
                unpaused.clear()

        fileExists = exists(outDir + "/" + fileName)
        outFile = open(outDir+ "/" + fileName, "a+")

        #write header if file did not exist yet
        if (not fileExists):
            outFile.write(header)
            outFile.write("\n")

        outFile.write(line)
        outFile.write("\n")

        outFile.close()

def getReference(deNovoDir):

    #get exons and splice junctions from de novo splice graphs
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    allGeneInfo = {}

    for gene in listdir(deNovoDir):

        #open file
        file = open(deNovoDir + "/" + gene, 'r')

        #get basic gene info from first line
        geneInfo = GeneInfo(file.readline())

        #Do not include unstranded genes
        if ((geneInfo.strand != "+") and (geneInfo.strand != "-")):
            continue

        while(True):
            line = file.readline()
            if not line: #end of file reached
                break

            #each line represents edge from one exon to another
            geneInfo.addJunction(line)

        #close file
        file.close()

        #add exons to "genome reference"
        for exon in geneInfo.observableExons:
            estart, estop = exon.split("-")

            #HTSeq uses 0-based positions and half-open interval (end position is not included)
            exons[HTSeq.GenomicInterval(geneInfo.chr, int(estart) - 1, int(estop))] += geneInfo.geneID

        #save gene info
        allGeneInfo[geneInfo.geneID] = geneInfo

    return exons, allGeneInfo

def readRawCounts(rawCountFile, allGeneInfo, bamfiles):

    #open file
    file = open(rawCountFile, 'r')

    #get cells from first line
    line = file.readline().rstrip("\n")
    cells = line.split('\t')

    #check if cells exist
    unalignedCells = []
    for cell in cells:
        if cell not in bamfiles:
            unalignedCells.append(cell)
    #error if none of the cells found
    if len(unalignedCells) == len(cells):
        raise RuntimeError("\n\nError: None of the cells in the count file had a corresponding bam file. \n"
                        "please check your if the header corresponds with the bam filenames (including .bam)")
    # #warning if few of cells missing
    # if len(unalignedCells) > 0:
    #     print("Warning: Some cells in the count file did not have a corresponding bam file, please check file names.")
    #     print("the following files are missing:")
    #     for cell in unalignedCells:
    #         print(cell)

    while (True):

        line = file.readline()
        if not line:
            break
        line = line.rstrip("\n")

        #read count line
        counts = line.split('\t')
        geneID = counts[0]

        #get gene information
        geneInfo = allGeneInfo.get(geneID)
        if geneInfo == None:
            continue

        #add length corrected read depth for each cell
        for i in range(1,len(counts)):
            if int(counts[i]) > 0:
                depth = float(counts[i])/int(geneInfo.longestTranscriptLen)
                geneInfo.addCellDepth(cells[i-1], depth)

def assignReadToGene(read, exons, cell, allGeneInfo):
    iset = None

    for cigop in read.cigar:
        #skip parts of alignment that are not alignment match, sequence match or sequence mismatch
        if not (cigop.type == "M" or cigop.type == "=" or cigop.type == "X") :
            continue
        for iv,step_set in exons[cigop.ref_iv].steps():
            if (len(step_set) > 0):
                if iset is None: #add to empty set
                    iset = step_set.copy()
                else:
                    iset.intersection_update(step_set) #only keep intersection of elements

    #read does not map to any gene
    if iset==None:
        return {}

    #get proportion to which read maps to each gene (proportional to raw gene counts)
    mappingGenes = {}
    sumGeneDepth = 0
    for gene in iset:
        depth = allGeneInfo.get(gene).getDepth(cell)
        if depth > 0:
            mappingGenes[gene] = depth
            sumGeneDepth += depth

    #normalize depths
    for gene in mappingGenes:
        mappingGenes[gene] /= sumGeneDepth

    return mappingGenes

def assignReadsToGene(read1, read2, exons, cell, allGeneInfo):

    iset = None

    #assign read 1
    for cigop in read1.cigar:
        #skip parts of alignment that are not alignment match, sequence match or sequence mismatch
        if not (cigop.type == "M" or cigop.type == "=" or cigop.type == "X") :
            continue
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
        for iv,step_set in exons[cigop.ref_iv].steps():
            if (len(step_set) > 0):
                if iset is None: #add to empty set
                    iset = step_set.copy()
                else:
                    iset.intersection_update(step_set) #only keep intersection of elements

    #read does not map to any gene
    if iset==None:
        return {}

    #get proportion to which read maps to each gene (proportional to raw gene counts)
    mappingGenes = {}
    sumGeneDepth = 0
    for gene in iset:
        depth = allGeneInfo.get(gene).getDepth(cell)
        if depth > 0:
            mappingGenes[gene] = depth
            sumGeneDepth += depth

    #normalize depths
    for gene in mappingGenes:
        mappingGenes[gene] /= sumGeneDepth

    return mappingGenes

def readPairedBam(sample, exons, allGeneInfo, bamDir, exonDepthQueue, junctionDepthQueue, geneCountQueue, lowQualQueue):

    #only start reading  bam when unpaused
    unpaused.wait()

    #gene specific coverage
    allGeneCov = {}

    #read in bam file to get aligned reads
    try:
        bam_reader = HTSeq.BAM_Reader(bamDir + "/" + sample)
        paired_bam_reader = HTSeq.pair_SAM_alignments_with_buffer(bam_reader)
    except Exception as e:
        print(bamDir + "/" + sample + " cannot be read")
        raise e

    for read1, read2 in paired_bam_reader:

        #skip unpaired reads & non primary alignments
        if (read1 is None or read2 is None or read1.not_primary_alignment or read2.not_primary_alignment):
            continue

        try:

            #get set of genes to which read maps
            mappingGenes = assignReadsToGene(read1, read2, exons, sample, allGeneInfo)

            #skip reads that don't map to any observed gene
            if (len(mappingGenes) == 0):
                continue

            #assign read (proportionally) to genes
            for gene in mappingGenes:

                #first time encountering gene?
                if gene not in allGeneCov:
                    allGeneCov[gene] = GeneCov(allGeneInfo.get(gene))

                depth = mappingGenes[gene]

                #add read to coverage vectors
                allGeneCov[gene].addPaired(read1, read2, depth)

        except KeyError:
            #chrom of read not in reference
            pass

    #write results to file
    for gene in allGeneCov:
        writeGeneDepthToFile(sample, allGeneInfo.get(gene), allGeneCov.get(gene), exonDepthQueue, junctionDepthQueue, geneCountQueue, lowQualQueue)

def readSingleBam(sample, exons, allGeneInfo, bamDir, exonDepthQueue, junctionDepthQueue, geneCountQueue, lowQualQueue):

    #only start reading  bam when unpaused
    unpaused.wait()

    #gene specific coverage
    allGeneCov = {}

    #read in bam file to get aligned reads
    try:
        bam_reader = HTSeq.BAM_Reader(bamDir + "/" + sample)

    except Exception as e:
        print(bamDir + "/" + sample + " cannot be read")
        raise e

    for i, read in enumerate(bam_reader):

        #skip non primary alignments
        if read.not_primary_alignment:
            continue

        try:

            #get set of genes to which read maps
            mappingGenes = assignReadToGene(read, exons, sample, allGeneInfo)

            #skip reads that don't map to any observed gene
            if (len(mappingGenes) == 0):
                continue

            #assign read (proportionally) to genes
            for gene in mappingGenes:

                #first time encountering gene?
                if gene not in allGeneCov:
                    allGeneCov[gene] = GeneCov(allGeneInfo.get(gene))

                depth = mappingGenes[gene]

                #add read to coverage vectors
                allGeneCov[gene].addSingle(read, depth)

        except KeyError:
            #chrom of read not in reference
            pass

    #write results to file
    for gene in allGeneCov:
        writeGeneDepthToFile(sample, allGeneInfo.get(gene), allGeneCov.get(gene), exonDepthQueue, junctionDepthQueue, geneCountQueue, lowQualQueue)

def writeGeneDepthToFile(sample, geneInfo, geneCov, exonDepthQueue, junctionDepthQueue, geneCountQueue, lowQualQueue):

    geneID = geneInfo.geneID

    #write exon counts to file
    line = sample
    header = ""
    for exon in geneInfo.getOrderedExons():
        #get genomic coordinates of exon
        estart, estop = exon.split("-")
        iv = HTSeq.GenomicInterval(geneInfo.chr, int(estart)-1, int(estop)) #HTSeq is 0-based and uses half-open interval

        sumExonDepth = sum(list(geneCov.exonCov[iv]))
        line += "\t" + str(sumExonDepth)
        header += "\t" + exon

    exonDepthQueue.put((geneID, header, line))

    #write junction counts to file
    line = sample
    header = ""

    for junction in geneInfo.getOrderedJunctions():
        #change from 0-based to 1-based and save intron positions instead of exon pos => +2 for start position
        start = junction[0] + 2     #change from 0-based to 1-based + save intron positions instead of exon pos
        stop = junction[1]          #change from 0-based to 1-based + save intron positions instead of exon pos

        depth = geneCov.junctionCounts.get(junction)
        line += "\t" + str(depth)
        header += "\t" + str(start) + "_" + str(stop)

    junctionDepthQueue.put((geneID, header, line))

    #write expressed genes to file
    geneCountQueue.put((geneID, "cell\tcounts", sample + "\t" + str(geneCov.highQualCounts)))

    #write proportion of low quality counts to file
    lowQualQueue.put((geneID,  "cell\tlowQualFrac", sample + "\t" + str(geneCov.getLowQualFrac())))

#remove non-expressed exons from depthfile
def filterDepth(depthFile):

    #read depth file
    depth = pd.read_table(depthFile, header = 0, index_col = 0)

    #filter exons/junctions not expressed in any cell
    depth = depth.loc[:,(depth>0).sum(axis=0)>0]

    #sort rows using alphabetical order of cellnames
    depth = depth.sort_index()

    #remove depthfile if no depth is observed
    if len(depth.columns) == 0:
        os.remove(depthFile)
    else:
        #replace file with filtered depth
        depth.to_csv(depthFile, sep = "\t")

def orderCells(file):
    #read file
    f = pd.read_table(file, header = 0, index_col = 0)

    #sort rows using alphabetical order of cellnames
    f = f.sort_index()

    #replace file with sorted file
    f.to_csv(file, sep = "\t")


def _parse_arguments():
    pa = argparse.ArgumentParser(
        add_help=False,
    )
    args, argv = pa.parse_known_args()

    pa = argparse.ArgumentParser(
        parents=[pa],
    )

    pa.add_argument(
        "--deNovo_dir",
        dest="deNovoDir",
        type=str,
        help="path/to/deNovoDir containing the splice graphs from the De Novo step",
    )

    pa.add_argument(
        "--bam_dir",
        dest="bamDir",
        type=str,
        help="path/to/bamDir containing the STAR output bam files",
    )

    pa.add_argument(
        "--geneCount_dir",
        dest="geneCountDir",
        type=str,
        help="path/to/geneCountDir where the gene counts will be written",
    )

    pa.add_argument(
        "--exonDepth_dir",
        dest="exonDepthDir",
        type=str,
        help="path/to/exonDepthDir where the output gene-specific exon depth files will be written"
    )

    pa.add_argument(
        "--junctionDepth_dir",
        dest="junctionDepthDir",
        type=str,
        help="path/to/junctionDepthDir where the output gene specific junction counts will be written"
    )

    pa.add_argument(
        "--lowQual_dir",
        dest="lowQualDir",
        type=str,
        help="path/to/lowQualDir where the fraction of low quality alignments will be writter"
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
        "--rawCountFile",
        dest="rawCountFile",
        type=str,
        help="count file containing raw gene counts (tab separated)"
    )

    pa.add_argument(
        "--selectCellsFile",
        dest="selectCellsFile",
        type=str,
        default = "",
        help="cells for which depth needs to be computed, if not provided all cells in bam_dir are used"
    )

    pa.add_argument(
        "--singleEnd",
        dest = "singleEnd",
        default=False,
        action='store_true',
        help="single end reads"
    )

    args = pa.parse_args()

    return args

def main():

    args = _parse_arguments()

    try:
        #read exons from de novo graphs as reference
        print("reading in exons and junctions as reference")
        exons, allGeneInfo = getReference(args.deNovoDir)

        #get all samples
        allBamfiles = [file for file in listdir(args.bamDir) if file.endswith(".bam")]

        #subset to cells to consider (if file provided)
        if args.selectCellsFile != "":
            bamfiles = []
            file = open(args.selectCellsFile, "r");
            bam = file.readline().strip()
            while (bam):
                #only select files that exist
                if bam in allBamfiles:
                    bamfiles.append(bam)
                else:
                    print("skipped cell " + bam + " : not present in " + args.bamDir)
                bam = file.readline().strip()
            file.close()
        else:
            #use all samples if no list is given
            bamfiles = allBamfiles

        #read raw genecount information
        print("reading raw counts")
        readRawCounts(args.rawCountFile, allGeneInfo, bamfiles)

        print("computing depth from bamfiles")
        manager = Manager()
        exonDepthQueue = manager.Queue()
        junctionDepthQueue = manager.Queue()
        geneCountQueue = manager.Queue()
        lowQualQueue = manager.Queue()

        #event to pause reading new bams if too many messages still need to be written to files
        global unpaused
        unpaused = manager.Event()
        unpaused.set() #no pause at startup

        nCPU = min(args.nCPU, len(bamfiles)) #avoid idle threads
        pool = Pool(nCPU+4) #4 additional processes for writers

        #first put listeners to work
        pool.apply_async(writeFiles, (args.junctionDepthDir, junctionDepthQueue, False, ))         #write junction depths
        pool.apply_async(writeFiles, (args.exonDepthDir, exonDepthQueue, True, ))                  #write exon depths
        pool.apply_async(writeFiles, (args.geneCountDir, geneCountQueue, False, ))                 #write gene counts
        pool.apply_async(writeFiles, (args.lowQualDir, lowQualQueue, False, ))                     #write low quality fraction

        #read bams
        if args.singleEnd:
            for _ in tqdm.tqdm(pool.imap_unordered(partial(readSingleBam, exons = exons, allGeneInfo = allGeneInfo,
                                                           bamDir = args.bamDir, exonDepthQueue = exonDepthQueue,
                                                           junctionDepthQueue = junctionDepthQueue,
                                                           geneCountQueue = geneCountQueue,
                                                           lowQualQueue = lowQualQueue), bamfiles), total=len(bamfiles)):
                pass
        else : #paired end
            for _ in tqdm.tqdm(pool.imap_unordered(partial(readPairedBam, exons = exons, allGeneInfo = allGeneInfo,
                                                           bamDir = args.bamDir, exonDepthQueue = exonDepthQueue,
                                                           junctionDepthQueue = junctionDepthQueue,
                                                           geneCountQueue = geneCountQueue,
                                                           lowQualQueue = lowQualQueue), bamfiles), total=len(bamfiles)):
                pass

        #stop writing processes
        exonDepthQueue.put("stop")
        junctionDepthQueue.put("stop")
        geneCountQueue.put("stop")
        lowQualQueue.put("stop")

        pool.close()
        pool.join()

        #remove non-expressed exons from depth files
        print("filtering non-expressed exons")
        pool = Pool(nCPU)
        exonDepthFiles = [(args.exonDepthDir + "/" + file) for file in listdir(args.exonDepthDir)]
        for _  in tqdm.tqdm(pool.imap_unordered(partial(filterDepth), exonDepthFiles), total=len(exonDepthFiles)):
            pass
        pool.close()
        pool.join()

        #remove non-expressed junctions from depth files
        print("filtering non-expressed junctions")
        pool = Pool(nCPU)
        junctionDepthFiles = [(args.junctionDepthDir + "/" + file) for file in listdir(args.junctionDepthDir)]
        for _ in tqdm.tqdm(pool.imap_unordered(partial(filterDepth), junctionDepthFiles), total=len(junctionDepthFiles)):
            pass
        pool.close()
        pool.join()

        #order geneCounts
        print("ordering cells in geneCounts")
        pool = Pool(nCPU)
        geneCountFiles = [(args.geneCountDir + "/" + file) for file in listdir(args.geneCountDir)]
        for _ in tqdm.tqdm(pool.imap_unordered(partial(orderCells), geneCountFiles), total=len(geneCountFiles)):
            pass
        pool.close()
        pool.join()

        #order geneCounts
        print("ordering cells in lowQual")
        pool = Pool(nCPU)
        lowQualFiles = [(args.lowQualDir + "/" + file) for file in listdir(args.lowQualDir)]
        for _ in tqdm.tqdm(pool.imap_unordered(partial(orderCells), lowQualFiles), total=len(lowQualFiles)):
            pass
        pool.close()
        pool.join()

    except:
        sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write(
            "  [Exception type: %s, raised in %s:%d]\n"
            % (
                sys.exc_info()[1].__class__.__name__,
                os.path.basename(traceback.extract_tb(sys.exc_info()[2])[-1][0]),
                traceback.extract_tb(sys.exc_info()[2])[-1][1],
            )
        )
        sys.exit(1)


if __name__ == "__main__":
    main()