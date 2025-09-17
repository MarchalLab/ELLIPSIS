#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 11:54:29 2025

@author: idlab297
"""

import HTSeq
import pandas as pd
import glob
import os
import argparse


if __name__ == "__main__":

    # Parse the command line arguments
    parser = argparse.ArgumentParser(description="get SJ files based on bam alignment")

    # Add command line arguments
    parser.add_argument('-b', '--bamDir', type=str, required=True, help="directory containing aligned bam files")
    parser.add_argument('-z', '--zUMIsDir', type=str, required=True, help="directory containing zUMIs output (mainly joint SJ.out.tab files)")

    args = parser.parse_args()
    
    #check if bamdir exists
    if not os.path.isdir(args.bamDir):
        raise Exception("bamDir ", args.bamDir, " does not exist")
    if not len(glob.glob(args.zUMIsDir + "/*.SJ.out.tab")) > 0:
        raise Exception(args.zUMIsDir + " does not contain SJ.out.tab file")
    
    #get observed junctions for all cells
    dfs = [pd.read_table(file, header = None) for file in glob.glob(args.zUMIsDir + "/*.SJ.out.tab")]
    novelJunctions = pd.concat(dfs, ignore_index = True)
    novelJunctions.columns = ["contig", "start", "end", "strand", "motif", "annot", "uniq", "multi", "overhang"]
    #novelJunctions = novelJunctions[novelJunctions.annot == 0]
    novelJunctions = novelJunctions[novelJunctions.uniq > 0]
    novelJunctions.strand = novelJunctions.strand.replace({0:'.', 1 : '+', 2 : '-'})
    
    novelJunctions = novelJunctions.groupby(["contig", "start", "end", "strand", "motif", "annot"]).agg({
        'uniq': 'sum',
        'multi' : 'sum',
        'overhang': 'max'
    }).reset_index()

    for file in glob.glob(args.bamDir + "/*.bam"):
        bam = HTSeq.BAM_Reader(file)
        
        cell = os.path.basename(file)
        cell = cell.replace(".bam", "")
        print(cell)
        
        thisJunctions = novelJunctions.copy(deep = True)
        thisJunctions.uniq = 0
        thisJunctions.multi = 0
        
        for read in bam:
            if not read.aligned or read.not_primary_alignment:
                continue
            
            for cigop in read.cigar:
                if cigop.type != "N":
                    continue
                
                thisJunctions.loc[((novelJunctions.contig == cigop.ref_iv.chrom) & 
                                               (novelJunctions.start == cigop.ref_iv.start + 1) & 
                                               (novelJunctions.end == cigop.ref_iv.end)), "uniq"] += 1
        
        thisJunctions = thisJunctions[thisJunctions.uniq > 0]
        thisJunctions.strand = thisJunctions.strand.replace({'.':0, '+':1, '-':2})
        thisJunctions.to_csv(args.bamDir +"/" + cell + "SJ.out.tab", 
                             sep = "\t", header = False, index = False)
