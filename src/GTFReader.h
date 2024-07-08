
#ifndef ELLIPSIS_GTFREADER_H
#define ELLIPSIS_GTFREADER_H

#include <string>
#include "spliceGraph.h"

using namespace std;

//Utility functions to read gtfs
class GTFReader{
private:

    string gtfFile;

    //gtf fields
    string geneIdField = "gene_id";
    string geneNameField = "gene_name";
    string transcriptIDField = "transcript_id";

    //current gene info
    string chrName, strand, geneID, geneName;
    size_t geneStart, geneEnd;

    //current transcript info
    string transcriptID;

    //current exon info
    size_t exonStart, exonEnd;

    /**
     * read gene info from gtf
     * @param line
     **/
    void readGeneLine(const string& line);

    /**
     * read transcript info from gtf
     * @param line
     */
    void readTranscriptLine(const string& line);

    /**
     * read exon info from gtf
     * @param line
     */
    void readExonLine(const string& line);

public:

    explicit GTFReader(string gtfFile) : gtfFile(move(gtfFile)){ }

    /**
     * build splice graphs from reference gtf file and write to directory
     * @param outDir directory where graphs will be written
     * @param genes genes for which graphs have to be written, if empty: all genes in the gtf file are written
     * @param geneInfoFile output file where gene information (chr, start, stop, strand) is written
     */
    void writeGraphs(const string& outDir, vector<string> genes, const string& geneInfoFile);

    /**
     * build splice graphs from reference gtf file and write to directory
     * @param outDir directory where graphs will be written
     */
    void writeGraphs(const string& outDir){
        vector<string> emptyVector;
        writeGraphs(outDir, emptyVector, "");
    }

    /**
     * get list of transcripts in which each node/edge occurs
     * @param graphDir directory containing graphs
     * @param nodeMapping
     * @param edgeMapping
     */
    void getTranscriptMapping(const string& graphDir, map<SpliceGraph*, map<SGNode*, vector<string>>>& nodeMapping,
                              map<SpliceGraph*, map<SGEdge*, vector<string>>>& edgeMapping);
};


#endif //ELLIPSIS_GTFREADER_H
