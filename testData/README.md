# ELLIPSIS tutorial on test dataset

To test if ELLIPSIS runs on your system, you can run ELLIPSIS on a small toy example provided in testData.
This dataset has 5 cells, and only contains reads for gene TSPAN6.

The required input consists of :
* STAR alignment results (incl. SJ.filtered.tab, which is the merged and filtered file for SJ.out.tab)
* the closest neighbors for each cell, here all 5 cells are each others closest neighbors
* raw gene count matrix

After installation of ELLIPSIS run:
```bash
conda activate ELLIPSIS_env

cd ELLIPSIS/testData

mkdir outELLIPSIS
../build/src/ELLIPSIS \
            -gtf TSPAN6.gtf \
            -alignmentDir alignment \
            -neighborFile neighbors.tsv \
            -readLength 75 \
            -outDir outELLIPSIS \
            -CPU 1 \
            -wSrcSink 1 \
            -wFlow 6 \
            -wClust 4 \
            -maxIter 100 \
            -minDepth 10 \
            -countFile counts/scCounts.tsv \
            -maxPaths inf \
            -novelJunctionFile alignment/SJ.filtered.tab \
            -verbose \
            > outELLIPSIS/ELLIPSIS.log
                 
mamba deactivate
```

This should generate the same output as in ELLIPSIS_results.
