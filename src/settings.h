
#ifndef ELLIPSIS_SETTINGS_H
#define ELLIPSIS_SETTINGS_H

#include <string>

#include "util.h"

class Settings {

private:
    // input file arguments
    std::string GTFFile;
    std::string alignmentDir;
    std::string clusterFilename;
    std::string neighborFile;
    std::string outDir;
    std::string runName;
    std::string trueGTF;
    std::string trueClusterFlow;
    std::string countFile;
    std::string novelJunctionFile;
    std::string selectCellsFile;

    size_t minJunctionCount;
    size_t minJunctionCells;
    size_t numCPU;
    size_t wObs;
    size_t wSrcSink;
    size_t wFlow;
    size_t wClust;
    size_t maxIter;
    size_t chunkSize;
    size_t maxPaths;

    double minDepth;
    double maxLowQual;

    bool pairedEnd;

    vector<std::string> selectedGenes;

    /**
     * Print usage to stdout
     */
    static void printUsage();

public:
    /**
     * Constructor that parses program arguments
     * @param argc Argument count
     * @param argv Actual arguments
     */
    Settings(int argc, char** argv);

    const std::string& getGTFFile() const;

    const std::string& getAlignmentDir() const;

    const std::string& getOutDir() const;

    std::string getAnnotGraphsDir() const;

    std::string getGeneCountDir() const;

    std::string getLowQualDir() const;

    std::string getGeneInfoFile() const;

    std::string getTrueGeneInfoFile() const;

    std::string getNovelSJDir() const;

    std::string getDeNovoDir() const;

    std::string getJunctionCountDir() const;

    std::string getSimplifiedDir() const;

    std::string getExonDepthDir() const;

    std::string getPSIDir() const;

    std::string getComparisonDir() const;

    std::string getClusterFile() const;

    std::string getNeighborFile() const;

    std::string getAlphaDir() const;

    std::string getAlphaLogDir() const;

    std::string getPSILogDir() const;

    std::string getFilteredDir() const;

    std::string getCountFile() const;

    std::string getNumNBDir() const;

    std::string getIterLog() const;

    std::string getComplexLog() const;

    const size_t& getMinJunctionCount() const;

    const size_t& getMinJunctionCells() const;

    const size_t& getNumCPU() const;

    const size_t& getWObs() const;

    const size_t& getWSrcSink() const;

    const size_t& getWFlow() const;

    const size_t& getWClust() const;

    const size_t& getMaxIter() const;

    const size_t& getMaxPaths() const;

    const size_t& getChunkSize() const;

    const double& getMinDepth() const;

    const double& getMaxLowQual() const;

    const string& getTrueGTF() const;

    string getTrueGraphsDir() const;

    string getTrueCellPSI() const;

    string getTrueClusterFlow() const;

    const vector<string>& getSelectedGenes() const;

    const string& getSelectCellsFile() const;

    const string& getNovelJunctionFile() const;

    bool isPairedEnd() const;
};


#endif //ELLIPSIS_SETTINGS_H
