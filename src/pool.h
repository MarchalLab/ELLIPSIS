
#ifndef ELLIPSIS_POOL_H
#define ELLIPSIS_POOL_H

#include <mutex>
#include <condition_variable>
#include <queue>
#include <map>

#include "util.h"
#include "SJCounter.h"

using namespace std;

//-----------------//
// Utility classes //
//-----------------//

class SJReader{
private:
    fstream ifs;

public:
    //current SJ
    string chrName, strand, junctionName;
    int start, end, annot, count;

    explicit SJReader(const string& SJFile);
    void close();
    bool readLine();
};

class SJCountReader{
private:
    fstream ifs;

public:
    //current SJ
    string junctionName;
    vector<int> counts;

    explicit SJCountReader(const string& SJCountFile);
    void close();
    bool readLine();

};

class PSIReader{
private:
    fstream ifs;

    void close();
    bool readLine();

public:
    vector<string> varNames;
    map<string, map<string, double>> PSI;

    explicit PSIReader(const string& fileName);
    map<string, map<string, double>> readAll();

};

class ThreadSafeStringSet{
private:
    mutex m;
    set<string> s;
public:
    ThreadSafeStringSet() = default;

    /**
     * add element to set (in thread safe manner)
     */
    void add(const string& to_add){
        lock_guard<mutex> lock(m);
        s.insert(to_add);
    }

    set<string> getSet(){
        lock_guard<mutex> lock(m);
        return s;
    }

    size_t getSize(){
        lock_guard<mutex> lock(m);
        return s.size();
    }
};

class Counter{
private:
    size_t ctr;
    mutex m;
public:
    Counter() : ctr(0){};

    /**
     * increment counter with 1
     */
    void increment(){
        lock_guard<mutex> lock(m);
        ctr += 1;
    }

    /**
     * increment counter with n
     * @param n
     */
    void increment(size_t n){
        lock_guard<mutex> lock(m);
        ctr += n;
    }

    /**
     * get counter value
     * @return
     */
    size_t getValue(){
        lock_guard<mutex> lock(m);
        return ctr;
    }
};

//------------------------------//
// Different types of arguments //
//------------------------------//

struct Args{
    Progressbar& progressbar;

    /**
     * Constructor
     * @param progressbar reference to progressbar
     */
    explicit Args(Progressbar& progressbar);

    /**
     * Do work for specific type of argument
     */
    virtual void doWork();

};


struct SJArgs : public Args{

    string SJFile;
    string cellID;
    map<string, SJCounter*> countMap;
    vector<GeneInfo> genes;
    ThreadSafeStringSet& junctionsWithoutGene;
    ThreadSafeStringSet& junctionsMultGenes;
    vector<string>& novelJunctions;

    /**
     * Constructor
     * @param progressbar reference to progressbar
     * @param cellID ID of cell corresponding to SJFile
     * @param countMap map containing gene and corresponding SJCounter
     * @param SJFile file to count SJ from
     * @param genes list with geneInfo for all considered genes
     * @param junctionsWithoutGene list of junctions that could not be attributed to any gene
     * @param junctionsMultGenes list of junctions that map could be attributed to multiple genes
     * @param novelJunctions novel junctions added in STAR pass2
     */
    SJArgs(Progressbar& progressbar, string cellID, map<string, SJCounter*>& countMap, string SJFile,
           vector<GeneInfo> genes, ThreadSafeStringSet& junctionsWithoutGene, ThreadSafeStringSet& junctionsMultGenes,
           vector<string>& novelJunctions);

    /**
     * aggregate SJ counts per gene
     */
    void doWork() override;
};

struct SJMergeArgs : public Args{

    string SJDir;

    Counter& numGenes;
    Counter& numJunctions;

    /**
     * Constructor
     * @param progressbar reference to progressbar
     * @param SJDir directory containing SJ counts per batch
     * @param numGenes count number of genes with novel SJs
     * @param numJunctions count number of novel SJs
     */
    SJMergeArgs(Progressbar& progressbar, string SJDir, Counter& numGenes, Counter& numJunctions);

    /**
     * merge partial SJ counts per gene
     */
    void doWork() override;
};

struct DeNovoArgs : public Args{
    string graphFile;
    string deNovoSJCountFile;
    string outFile;
    size_t minJunctionCount;
    size_t minCellCount;
    Counter& numDeNovo;
    Counter& numGenesChanged;

    /**
     * constructor
     * @param progressbar reference to progressbar
     * @param graphFile original splice graph (from reference)
     * @param deNovoSJCountFile SJcount file (this file might not exist if no novel junctions are observed for this gene)
     * @param outFile file where expanded splice graph will be written
     * @param minJunctionCount minimum number of counts for a new junction to be accepted
     * @param minCellCount minimum number of cells for a new junction occurs to be accepted
     * @param numDeNovo counter keeping track of De Novo junctions added
     * @param numGenesChanged counter keeping track of number of genes for which de novo junction is added
     */
    DeNovoArgs(Progressbar & progressbar, string graphFile, string deNovoSJCountFile, string outFile, size_t minJunctionCount,
               size_t minCellCount, Counter& numDeNovo, Counter& numGenesChanged);

    /**
     * Add novel observed junctions to the splice graphs
     */
    void doWork() override;
};

struct SimplifyArgs : public Args{

    string graphFile;
    string exonDepthFile;
    string junctionCountFile;
    string geneCountFile;
    string outFile;

    double minDepth;
    size_t minJunctionCount;
    size_t minJunctionCells;

    Counter& numNodesRemoved;
    Counter& numGraphsNodesRemoved;
    Counter& numEdgesRemoved;
    Counter& numGraphsEdgesRemoved;
    Counter& numMerged;
    Counter& numGraphsMerged;
    Counter& numNodesRemaining;
    Counter& numEdgesRemaining;

    /**
     * Constructor
     * @param progressbar reference to progressbar
     * @param graphFile file containing splice graph to be simplified
     * @param exonDepthFile file containing depth per observed exon and per cell
     * @param junctionCountFile file containing junction counts per observed junction and per cell
     * @param geneCountFile file containing gene counts per cell
     * @param minDepth minimum number depth for a gene to be considered, genes that are lowly expressed in all cells are discarded here
     * @param minJunctionCount minimum depth for an observed junction with unobserved start/end exon(s) to be considered expressed
     * @param minJunctionCells minimum number of cells where an observed junction with unobserved start/end exon(s) is expressed to be included in the splicegraph
     * @param outFile output file to which simplified graph is written
     * @param numNodesRemoved keep track of number of nodes removed
     * @param numGraphsNodesRemoved keep track of number of graphs where at least 1 node is removed
     * @param numEdgesRemoved keep track of number of edges removed
     * @param numGraphsEdgesRemoved keep track of number of graphs where at least 1 edge is removed
     * @param numMerged keep track of number of merge operations performed
     * @param numGraphsMerged keep track of number of graphs where at least 1 merge operation is performed
     * @param numNodesRemaining keep track of number of nodes in the simplified graph
     * @param numEdgesRemaining keep track of number of edges in the simplified graph
     */
    SimplifyArgs(Progressbar& progressbar, string graphFile, string exonDepthFile, string junctionCountFile,
                 string geneCountFile, double minDepth, size_t minJunctionCount, size_t minJunctionCells,
                 string outFile, Counter& numNodesRemoved, Counter& numGraphsNodesRemoved,
                 Counter& numEdgesRemoved, Counter& numGraphsEdgesRemoved,  Counter& numMerged,
                 Counter& numGraphsMerged, Counter& numNodesRemaining, Counter& numEdgesRemaining);


    /**
     * simplify graphs by only using observed exons & junctions (while keeping it connected)
     */
    void doWork() override;
};

struct PSIArgs : public Args{
    string graphFile;
    string exonDepthFile;
    string junctionCountFile;
    string outFile;
    string alphaFile;
    string alphaLogDir;
    string PSILogDir;
    string filteredFile;
    string lowQualFile;
    string numNBFile;

    map<string, vector<string>>& neighbors;

    size_t wObs;
    size_t wSrcSink;
    size_t wFlow;
    size_t wClust;
    size_t maxIter;
    uint1024_t maxPaths;

    double minDepth;
    double maxLowQual;

    LogWriter& nIterLog;
    LogWriter& complexLog;

    /**
     * Constructor
     * @param progressbar reference to progressbar
     * @param nIterLog logging the nr of iterations needed until convergence per file
     * @param complexLog logging the genes that are ignored because they are too complex
     * @param graphFile file containing splice graph for which PSI values are computed
     * @param exonDepthFile file containing exon depths for this gene
     * @param junctionCountFile file containing junction counts for this gene
     * @param outFile file to write computed PSI values to
     * @param alphaFile  file to write computed gene depth (alpha) to
     * @param alphaLogDir directory to write intermediate alphas to (only used in verbose mode)
     * @param PSILogDir directory to write intermediate PSI results to (only used in verbose mode)
     * @param filteredFile file to write filtered cells information to
     * @param lowQualFile file containing low quality read fractions for each cell
     * @param numNBFile file containing number of (expressed) neighbors for each cell
     * @param neighbors list of n neighboring cells
     * @param wObs weight assigned to equations concerning the observed depth
     * @param wSrcSink weight assigned to equations that impose PSI in src/sink = 100%
     * @param wFlow weight assigned to equations concerning the conservation of flow of PSI values
     * @param wClust weight assigned to equations concerning intra-cell type similarity
     * @param minDepth minimum depth (initial value for alpha) needed for cell to compute PSI values
     * @param maxIter maximum number of iterations in PSI computation
     * @param maxLowQual maximum fraction of low quality reads for cell to be considered
     * @param maxPaths maximum number of source to sink paths to compute PSI values
     */
    PSIArgs(Progressbar &progressbar, LogWriter &nIterLog, LogWriter &complexLog, string graphFile,
            string exonDepthFile, string junctionCountFile, string outFile, string alphaFile, string alphaLogDir,
            string PSILogDir, string filteredFile, string lowQualFile, string numNBFile,
            map<string,vector<string>>& neighbors, size_t wObs, size_t wSrcSink, size_t wFlow, size_t wClust, double minDepth,
            size_t maxIter, double maxLowQual, uint1024_t maxPaths);

    /**
     * compute PSI values for each exon/junction in this gene
     */
    void doWork() override;
};

struct CompareArgs : public Args {

    string trueGraphFile;
    string truePSIFile;
    string estGraphFile;
    string estPSIFile;

    string outDir;

    /**
     * Constructor
     * @param progressbar reference to progressbar
     * @param trueGraphFile splice graph for ground truth
     * @param truePSIFile PSI values for ground truth
     * @param estGraphFile splice graph for estimated PSI values
     * @param estPSIFile estimated PSI values
     * @param outDir output directory to write files to
     */
    CompareArgs(Progressbar &progressbar, string trueGraphFile, string truePSIFile, string estGraphFile,
                string estPSIFile, string outDir);

    /**
     * Compute accuracy scores for one gene/graph
     */
    void doWork() override;
};

//------------------------------//
// implementation of threadpool //
//------------------------------//

class Pool {
private:
    mutex pool_mutex;
    condition_variable condition;
    queue<Args*> jobQueue;
    bool finished;

public:
    /**
     * Default constructor
     */
    Pool();

    /**
     * Add new job to pool
     * @param arg arguments for mergeCounts function
     */
    void addJob(Args* args);

    /**
     * Infinite loop checking for new jobs, and executing them
     */
    void doWork();

    /**
     * Terminate pool:
     */
    void terminatePool();
};


#endif //ELLIPSIS_POOL_H



