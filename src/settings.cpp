
#include <iostream>
//#include <sys/stat.h>

#include "settings.h"
#include "util.h"

using namespace std;

//global variables
size_t READLENGTH;
bool VERBOSE = false;

void Settings::printUsage()
{
    cout <<
        "Usage: ELLIPSIS -gtf ref.gtf -alignmentDir /path/to/alignment -clusterFile clusters.tsv -readLength <RL> [options]\n\n"

        " required arguments\n"
        "  -gtf         \t\tref.gtf the gtf file (from Ensembl) containing the reference genome annotation\n"
        "  -alignmentDir\t\t/path/to/alignment the directory containing the STAR output files (bam files and SJ.out.tab files)\n"
        "  -neighborFile\t\tneighbors.tsv the file containing a list of k nearest neighbors (tab separated)\n"
        "  -countFile   \t\tcounts.tsv the file containing the geneCounts (tab separated)\n"
        "  -readLength  \t\tthe average length of a read (for PE reads: length of 1 mate)\n\n"

        " [addtional options without argument]\n"
        "  -singleEnd  \t\tsingle end reads [default PE reads]"
        "  -verbose    \t\tcreate additional log files during PSI calculations (logPSI and logAlpha) \n\n"

        " [additional options with 1 argument]\n"
        "  -novelJunctionFile\t\tnovel junctions added in STAR pass2\n"
        "  -minJunctionCount \t\tminimum number of junction spanning reads per cell required for a non annotated junction to be accepted [default = 5]\n"
        "  -minJunctionCells \t\tminimum number of cells for which non annotated junction is expressed for it to be accepted [default = 10]\n"
        "  -minDepth         \t\tminimum read depth required to consider a cell for computing PSI values [default = 10]\n"
        "  -wObs             \t\tweight assigned to observed depths [default = 1]\n"
        "  -wSrcSink         \t\tweight assigned to impose src/sink have PSI = 100% (high values are advised to avoid PSI-values > 100%) [default = 1]\n"
        "  -wFlow            \t\tweight assigned to impose conservation of flow of PSI values [default = 6]\n"
        "  -wClust           \t\tweight assigned to favor intra-cell type similarity [default = 4]\n"
        "  -maxIter          \t\tmaximum number of iterations in compute PSI step [default = 100]\n"
        "  -maxLowQual       \t\tmaximum fraction of low quality read mappings to gene [default = 0.1]\n"
        "  -maxPaths         \t\tmaximum number of source to sink paths to compute PSI values (use inf to remove maxPaths filter) [default = 1e9]\n"
        "  -chunkSize        \t\tmaximum number of files read at once for SJMerge [default = 1000]\n"

        "  -outDir           \t\toutput directory [default = this directory]\n"
        "  -run              \t\trun name \n"
        "  -CPU              \t\tnumber of CPU threads [default = 4]\n"
        "  -selectGenes      \t\tfile containing genes to consider. If not provided, flow is computed for all genes possible\n"
        "  -selectCells      \t\tfile containing subset of cells to consider. If not provided, all cells are considered\n"


        "  -trueClusterFlow  \t\tdirectory containing trueClusterFlow files \n"
        "  -clusterFile      \t\tclusters.tsv the file containing a list of cells and their corresponding cluster (tab separated)\n"
        "  -trueGTF          \t\tgtf file containing only transcripts that are present, if option not given, no trueGraphs are build.\n\n";

}

Settings::Settings(int argc, char ** argv) : minJunctionCount(5), outDir(), runName(),
    numCPU(4), wObs(1), wSrcSink(1), wFlow(6), wClust(4), maxIter(100), minDepth(10), minJunctionCells(10),
    trueGTF(), trueClusterFlow(), maxLowQual(0.1), novelJunctionFile(), maxPaths(1e9), chunkSize(1000), pairedEnd(true)
{
    const size_t reqArguments = 4;

    // no arguments provided
    if (argc < reqArguments) {
        printUsage();
        exit(EXIT_FAILURE);
    }

    for (size_t i = 1; i < argc; i++) {
        string arg(argv[i]);

        // process options without arguments
        if (arg == "-verbose"){
            VERBOSE = true;
            continue; //no arg
        }
        if (arg == "-singleEnd"){
            pairedEnd = false;
            continue; //no arg
        }

        try {
            // process options with 1 argument
            if (arg == "-gtf") {
                GTFFile = argv[i + 1];
            } else if (arg == "-alignmentDir") {
                alignmentDir = argv[i + 1];
            } else if (arg == "-minJunctionCount") {
                minJunctionCount = Util::string2uint(argv[i + 1]);
            } else if (arg == "-minJunctionCells") {
                minJunctionCells = Util::string2uint(argv[i + 1]);
            } else if (arg == "-outDir") {
                outDir = argv[i + 1];
                if (!Util::dirExists(outDir))
                    Util::createDir(outDir);
            } else if (arg=="-run"){
                runName = argv[i+1];
            } else if (arg == "-CPU") {
                numCPU = Util::string2uint(argv[i + 1]);
            } else if (arg == "-wObs") {
                wObs = Util::string2uint(argv[i + 1]);
                if (wObs == 0)
                    throw invalid_argument("value cannot be zero");
            } else if (arg == "-wSrcSink") {
                wSrcSink = Util::string2uint(argv[i + 1]);
            } else if (arg == "-wFlow") {
                wFlow = Util::string2uint(argv[i + 1]);
                if (wFlow == 0)
                    throw invalid_argument("value cannot be zero");
            } else if (arg == "-wClust") {
                wClust = Util::string2uint(argv[i + 1]);
            } else if (arg == "-maxIter") {
                maxIter = Util::string2uint(argv[i + 1]);
            } else if (arg == "-minDepth") {
                minDepth = stod(argv[i + 1]);
            } else if (arg == "-maxLowQual") {
                maxLowQual = stod(argv[i + 1]);
            } else if (arg == "-trueGTF") {
                trueGTF = argv[i + 1];
            } else if (arg == "-trueClusterFlow") {
                trueClusterFlow = argv[i + 1];
            } else if (arg == "-clusterFile") {
                clusterFilename = argv[i + 1];
            } else if (arg == "-readLength") {
                READLENGTH = Util::string2uint(argv[i + 1]);
                if (READLENGTH == 0)
                    throw invalid_argument("value cannot be zero");
            } else if (arg == "-selectGenes") {
                string selectGeneFile = argv[i + 1];
                if (!Util::fileExists(selectGeneFile))
                    throw runtime_error("please provide a valid list of genes, " + selectGeneFile + " does not exist");
                fstream genesIFS(selectGeneFile.c_str());
                string line;
                while (getline(genesIFS, line)) {
                    selectedGenes.push_back(line);
                }
            } else if (arg == "-selectCells") {
                selectCellsFile = argv[i + 1];
                if (!Util::fileExists(selectCellsFile))
                    throw runtime_error("please provide a valid list of cells, " + selectCellsFile + " does not exist");
            } else if (arg == "-neighborFile") {
                neighborFile = argv[i + 1];
            } else if (arg == "-countFile") {
                countFile = argv[i + 1];
            } else if (arg == "-novelJunctionFile") {
                novelJunctionFile = argv[i + 1];
            } else if (arg == "-maxPaths") {
                if (Util::compareCaseInsensitiveStrings(argv[i + 1], "inf"))
                    maxPaths = 0; //maxpaths = 0 => no limit to max nr of paths
                else{
                    maxPaths = Util::string2uint(argv[i + 1]);
                    if (maxPaths == 0)
                        throw invalid_argument("value cannot be zero");
                }
            } else if (arg == "-chunkSize") {
                chunkSize = Util::string2uint(argv[i + 1]);
                if (chunkSize == 0)
                    throw invalid_argument("value cannot be zero");
            } else {
                cerr << "unknown argument " << arg << endl;
            }
        }catch(const invalid_argument& e){
            cerr << "please give valid value for " << argv[i] << endl;
            cerr << e.what() << endl;
            exit(EXIT_FAILURE);
        }

        i++; //skip argument
    }

    //double readlength if paired end reads
    if (pairedEnd)
        READLENGTH *= 2;

    //check if all required arguments are given
    if (GTFFile.empty()) {
        cerr << "Specify annotation file" << endl;
        printUsage();
        exit(EXIT_FAILURE);
    }
    if (alignmentDir.empty()){
        cerr << "Specify the path to the SJ.out.tab files" << endl;
        printUsage();
        exit(EXIT_FAILURE);
    }
    if (neighborFile.empty()) {
        cerr << "Specify neighbor file" << endl;
        printUsage();
        exit(EXIT_FAILURE);
    }if (countFile.empty()){
        cerr << "specify count file" << endl;
        printUsage();
        exit(EXIT_FAILURE);
    }

    //check if directories exist
    if (!Util::fileExists(GTFFile))
        throw runtime_error("cannot open file " + GTFFile);
    if (!Util::dirExists(alignmentDir)){
        throw runtime_error(alignmentDir + " does not exist");
    }
    if (!clusterFilename.empty() and !Util::fileExists(clusterFilename) )
        throw runtime_error(clusterFilename + " does not exist");
    if (!neighborFile.empty() and !Util::fileExists(neighborFile))
        throw runtime_error(neighborFile + " does not exist");
    if (!Util::fileExists(countFile))
        throw runtime_error(countFile + " does not exist");

    //check if all required files are provided for comparison with ground truth
    if ((! trueClusterFlow.empty()) | (! trueGTF.empty())){
        if ((trueClusterFlow.empty()) | (trueGTF.empty()) | (clusterFilename.empty()))
            throw runtime_error("please provide the true gtf file, the true cluster flow and the cell clustering information if you want to compare with the ground truth");
    }

    //create run specific directory
    if (! runName.empty())
        Util::createDir(outDir + "/" + runName);

}

const std::string& Settings::getGTFFile() const {
    return GTFFile;
}

const std::string& Settings::getAlignmentDir() const {
    return alignmentDir;
}

const std::string& Settings::getOutDir() const{
    return outDir;
}

std::string Settings::getAnnotGraphsDir() const{
    return outDir + "/annotGraphs";
}

std::string Settings::getGeneInfoFile() const{
    return outDir + "/geneInfo.tsv";
}

std::string Settings::getTrueGeneInfoFile() const{
    return outDir + "/trueGeneInfo.tsv";
}

std::string Settings::getNovelSJDir() const{
    return outDir + "/novelSJCounts";
}

std::string Settings::getDeNovoDir() const {
    return outDir + "/deNovoGraphs";
}

std::string Settings::getJunctionCountDir() const {
    return outDir + "/allJunctionCounts";
}

std::string Settings::getExonDepthDir() const {
    return outDir + "/exonDepth";
}

std::string Settings::getGeneCountDir() const{
    return outDir + "/geneCounts";
}

std::string Settings::getLowQualDir() const{
    return outDir + "/lowQual";
}

std::string Settings::getSimplifiedDir() const {
    return outDir + "/simplifiedGraphs";
}

std::string Settings::getPSIDir() const {
    if (runName.empty())
        return outDir + "/PSI";
    return outDir + "/" + runName + "/PSI";
}

std::string Settings::getComparisonDir() const {
    if (runName.empty())
        return outDir + "/compareGraphs";
    return outDir + "/" + runName + "/compareGraphs";
}

std::string Settings::getClusterFile() const {
    return clusterFilename;
}

std::string Settings::getNeighborFile() const {
    return neighborFile;
}

std::string Settings::getAlphaDir() const{
    if (runName.empty())
        return outDir + "/alpha";
    return outDir + "/" + runName + "/alpha";
}

std::string Settings::getAlphaLogDir() const{
    if (runName.empty())
        return outDir + "/logAlpha";
    return outDir + "/" + runName + "/logAlpha";
}

std::string Settings::getPSILogDir() const{
    if (runName.empty())
        return outDir + "/logPSI";
    return outDir + "/" + runName + "/logPSI";
}

std::string Settings::getFilteredDir() const {
    if (runName.empty())
        return outDir + "/filteredCells";
    return outDir + "/" + runName + "/filteredCells";
}

std::string Settings::getCountFile() const {
    return countFile;
}

std::string Settings::getNumNBDir() const{
    if (runName.empty())
        return outDir + "/numNeighbors";
    return outDir + "/" + runName + "/numNeighbors";
}


std::string Settings::getIterLog() const{
    if (runName.empty())
        return outDir + "/nIterPSI.log";
    return outDir + "/" + runName + "/nIterPSI.log";
}

std::string Settings::getComplexLog() const{
    if (runName.empty())
        return outDir + "/genesTooComplex.log";
    return outDir + "/" + runName + "/genesTooComplex.log";
}

const size_t& Settings::getMinJunctionCount() const {
    return minJunctionCount;
}

const size_t& Settings::getMinJunctionCells() const {
    return minJunctionCells;
}

const size_t& Settings::getNumCPU() const {
    return numCPU;
}

const size_t& Settings::getWObs() const{
    return wObs;
}

const size_t& Settings::getWSrcSink() const{
    return wSrcSink;
}

const size_t& Settings::getWFlow() const{
    return wFlow;
}

const size_t& Settings::getWClust() const{
    return wClust;
}

const size_t& Settings::getMaxIter() const {
    return maxIter;
}

const size_t& Settings::getMaxPaths() const {
    return maxPaths;
}

const size_t& Settings::getChunkSize() const {
    return chunkSize;
}

const double& Settings::getMinDepth() const {
    return minDepth;
}

const double& Settings::getMaxLowQual() const{
    return maxLowQual;
}

const string& Settings::getTrueGTF() const {
    return trueGTF;
}

string Settings::getTrueGraphsDir() const{
    return outDir + "/trueGraphs";
}

string Settings::getTrueCellPSI() const{
    return outDir + "/trueCellPSI";
}

string Settings::getTrueClusterFlow() const{
    return trueClusterFlow;
}

const vector<string>& Settings::getSelectedGenes() const{
    return selectedGenes;
}

const string& Settings::getSelectCellsFile() const{
    return selectCellsFile;
}

const string& Settings::getNovelJunctionFile() const{
    return novelJunctionFile;
}

bool Settings::isPairedEnd() const{
    return pairedEnd;
}