
#include "pool.h"
#include "util.h"
#include "observedDepth.h"
#include "PSIsolver.h"
//#include <sstream>
#include <iostream>
#include <cassert>

//global variables
extern int READLENGTH;
extern bool VERBOSE;

//------------------------------//
// Different types of arguments //
//------------------------------//

Args::Args(Progressbar& progressbar) : progressbar(progressbar) {}

void Args::doWork(){
    throw runtime_error("No work can be done for general Args. You need a specific type of arguments.");
}

SJArgs::SJArgs(Progressbar &progressbar, string cellID, map<string, SJCounter*>& countMap,
               string SJFile, vector<GeneInfo> genes, ThreadSafeStringSet& junctionsWithoutGene,
               ThreadSafeStringSet& junctionsMultGenes, vector<string>& novelJunctions) : Args(progressbar),
    cellID(move(cellID)), countMap(countMap), SJFile(move(SJFile)), genes(move(genes)), junctionsWithoutGene(junctionsWithoutGene),
    junctionsMultGenes(junctionsMultGenes), novelJunctions(novelJunctions){}

void SJArgs::doWork() {

    SJReader fileReader(SJFile);

    //get junction counts
    while (fileReader.readLine()){

        //do not consider junctions without any uniquely mapping reads
        if (fileReader.count == 0)
            continue;

        //do not consider already annotated junctions (that are not added in STAR pass2)
        if ((fileReader.annot == 1) &
            (std::find(novelJunctions.begin(), novelJunctions.end(), fileReader.junctionName) == novelJunctions.end()))
                continue;

        //find gene(s) to which junction belongs
        vector<string> geneAssignments;
        for (auto const &gene : genes) {
            if (gene.containsJunction(fileReader.chrName, fileReader.start, fileReader.end)){
                geneAssignments.push_back(gene.getGeneID());
            }
        }
        if (geneAssignments.empty()){ //junction cannot be attributed to any gene
            junctionsWithoutGene.add(fileReader.junctionName);
        }else{
            if (geneAssignments.size() >= 2) { //junction cannot be attributed to unique gene -> add to all of them
                junctionsMultGenes.add(fileReader.junctionName);
            }
            for (auto const &gene : geneAssignments){
                countMap[gene]->addCount(fileReader.junctionName, cellID, fileReader.count);
            }
        }

    }

    fileReader.close();

    progressbar.update(1);
}

SJMergeArgs::SJMergeArgs(Progressbar& progressbar, string SJDir, Counter& numGenes, Counter& numJunctions) :
    Args(progressbar), SJDir(move(SJDir)), numGenes(numGenes), numJunctions(numJunctions) {}

void SJMergeArgs::doWork(){

    //read counts in all batches
    SJCounter mergedCounts;
    for (const auto& batch : Util::getFileList(SJDir)){
        SJCounter newCounts(SJDir + "/" + batch);
        mergedCounts.merge(newCounts);
    }

    if (mergedCounts.getNumJunctions() > 0){
        numGenes.increment();
        numJunctions.increment(mergedCounts.getNumJunctions());
    }

    //remove directory with partial files
    Util::removeDir(SJDir);

    //write counts to file
    mergedCounts.writeToFile(SJDir);

    progressbar.update(1);
}

DeNovoArgs::DeNovoArgs(Progressbar & progressbar, string graphFile, string deNovoSJCountFile, string outFile, size_t minJunctionCount,
                       size_t minCellCount, Counter& numDeNovo, Counter& numGenesChanged) : Args(progressbar),
       graphFile(move(graphFile)), deNovoSJCountFile(move(deNovoSJCountFile)), outFile(move(outFile)),
       minJunctionCount(minJunctionCount), minCellCount(minCellCount), numDeNovo(numDeNovo), numGenesChanged(numGenesChanged) {}

void DeNovoArgs::doWork() {

    SpliceGraph graph(graphFile);
    bool graphChanged = false;

    //only do something if novel splice junctions are observed
    if (Util::fileExists(deNovoSJCountFile)){
        SJCountReader countReader(deNovoSJCountFile);

        while (countReader.readLine()){
            vector<int> counts = countReader.counts;

            //get nr of cells that express junction (at least minJunctionCount times)
            int nCellsExpr = 0;
            for (auto const &count : counts){
                if (count > minJunctionCount){
                    nCellsExpr ++;
                }
            }

            if (nCellsExpr > minCellCount){
                bool junctionAdded = graph.addJunction(countReader.junctionName);
                if (junctionAdded){
                    numDeNovo.increment();
                    graphChanged = true;
                }
            }
        }

        countReader.close();
    }

    if (graphChanged)
        numGenesChanged.increment();

    //write resulting graph to file
    ofstream countsOFS(outFile);
    countsOFS << graph;
    countsOFS.close();

    //remove all nodes and edges from graph
    graph.deleteGraph();

    progressbar.update(1); //update progressbar
}


SimplifyArgs::SimplifyArgs(Progressbar &progressbar, string graphFile, string exonDepthFile, string junctionCountFile,
                           string geneCountFile, double minDepth, size_t minJunctionCount, size_t minJunctionCells,
                           string outFile, Counter& numNodesRemoved, Counter& numGraphsNodesRemoved,
                           Counter& numEdgesRemoved, Counter& numGraphsEdgesRemoved, Counter& numMerged,
                           Counter& numGraphsMerged, Counter& numNodesRemaining,
                           Counter& numEdgesRemaining) : Args(progressbar),
        graphFile(move(graphFile)), exonDepthFile(move(exonDepthFile)), junctionCountFile(move(junctionCountFile)),
        geneCountFile(move(geneCountFile)), minDepth(minDepth), minJunctionCount(minJunctionCount),
        minJunctionCells(minJunctionCells), outFile(move(outFile)), numNodesRemoved(numNodesRemoved),
        numGraphsNodesRemoved(numGraphsNodesRemoved), numEdgesRemoved(numEdgesRemoved),
        numGraphsEdgesRemoved(numGraphsEdgesRemoved), numMerged(numMerged), numGraphsMerged(numGraphsMerged),
        numNodesRemaining(numNodesRemaining), numEdgesRemaining(numEdgesRemaining) {}

void SimplifyArgs::doWork() {

    //get graph
    SpliceGraph graph(graphFile);

    int minGeneExpr = floor(minDepth * graph.getShortestTranscriptLen()/READLENGTH);

    //check if gene is expressed enough
    bool geneExpressed = false;
    fstream countIFS(geneCountFile);
    string line, cellID; double cellCount;
    getline(countIFS, line); //skip first line (header)
    while (getline(countIFS, line)){
        istringstream iss(line);
        iss >> cellID >> cellCount;
        if (cellCount > minGeneExpr){
            geneExpressed = true;
            break;
        }
    }
    countIFS.close();

    if (! geneExpressed){
        graph.deleteGraph();
        progressbar.update(1);
        return;
    }

    //get observed exons
    set<SGNode*> observedNodes = {graph.getSourceNode(), graph.getSinkNode()}; //already add artificial nodes
    graph.getObservedNodes(exonDepthFile, observedNodes);

    //get observed edges
    set<SGEdge*> observedEdges;
    graph.getObservedEdges(junctionCountFile, observedEdges, observedNodes, minJunctionCount, minJunctionCells);

    //make sure graph of observednodes and observed edges is connected
    graph.makeConnected(observedNodes, observedEdges);

    //remove unobserved edges from graph
    bool edgesRemoved = false;
    vector<SGEdge*> edges = graph.getEdges();
    for (auto const &edge : edges){
        if (observedEdges.find(edge) == observedEdges.end()) {
            graph.removeEdge(edge);
            numEdgesRemoved.increment();
            edgesRemoved = true;
        }
    }

    if (edgesRemoved)
        numGraphsEdgesRemoved.increment();

    //remove unobserved nodes from graph
    bool nodesRemoved = false;
    vector<SGNode*> nodes = graph.getNodes();
    for (auto const &node : nodes){
        if (observedNodes.find(node) == observedNodes.end()){
            graph.removeNode(node);
            numNodesRemoved.increment();
            nodesRemoved = true;
        }
    }

    if(nodesRemoved)
        numGraphsNodesRemoved.increment();

    //merge consecutive nodes if possible
    size_t nMerged = graph.mergeNodes();
    numMerged.increment(nMerged);
    if (nMerged > 0)
        numGraphsMerged.increment();

    if (graph.getNumEdges() != 0 && graph.isSpliced()){ //only consider non-empty graphs, and genes that where splicing is possible

        //write resulting graph to file
        ofstream OFS(outFile);
        OFS << graph;
        OFS.close();

        //get number of remaining nodes and edges
        size_t numNodes = graph.getNumNodes();
        size_t numEdges = graph.getNumEdges();

        numNodesRemaining.increment(numNodes);
        numEdgesRemaining.increment(numEdges);
    }

    //reset graph
    graph.deleteGraph();

    //update progressbar
    progressbar.update(1);
}

PSIArgs::PSIArgs(Progressbar &progressbar, LogWriter &nIterLog, LogWriter &complexLog, string graphFile, string exonDepthFile,
                 string junctionCountFile, string outFile, string alphaFile, string alphaLogDir, string PSILogDir,
                 string filteredFile, string lowQualFile, string numNBFile, map<string, vector<string>>& neighbors,
                 size_t wObs, size_t wSrcSink, size_t wFlow, size_t wClust, double minDepth, size_t maxIter,
                 double maxLowQual, uint1024_t maxPaths) : Args(progressbar),
         nIterLog(nIterLog), complexLog(complexLog), graphFile(move(graphFile)), exonDepthFile(move(exonDepthFile)),
         junctionCountFile(move(junctionCountFile)), outFile(move(outFile)), alphaFile(move(alphaFile)),
         alphaLogDir(move(alphaLogDir)), PSILogDir(move(PSILogDir)), filteredFile(move(filteredFile)),
         lowQualFile(move(lowQualFile)), numNBFile(move(numNBFile)), neighbors(neighbors), wObs(wObs),
         wSrcSink(wSrcSink), wFlow(wFlow), wClust(wClust), minDepth(minDepth), maxIter(maxIter), maxLowQual(maxLowQual),
         maxPaths(move(maxPaths)) {}

void PSIArgs::doWork(){
    //get graph
    SpliceGraph* spliceGraph = new SpliceGraph(graphFile);
    VarGraph* graph = new VarGraph(spliceGraph);

    if (maxPaths > 0){
        //extremely complex graphs are ignored
        uint1024_t numPaths = graph->getTotalNumberOfPaths();
        if (numPaths > maxPaths) {
            std::stringstream numPathsString;
            numPathsString << numPaths;
            complexLog.writeLine(graph->getGeneID() + "\t" + numPathsString.str());

            //delete graphs
            spliceGraph->deleteGraph();
            delete (spliceGraph);
            graph->deleteGraph();
            delete (graph);

            //update progressbar
            progressbar.update(1);

            return;
        }

    }

    if (VERBOSE) {
        Util::createDir(PSILogDir);
        Util::createDir(alphaLogDir);
    }

    map<string, double> lowQualFrac = Util::readLowQualFile(lowQualFile);

    //get observed depth
    ObservedDepth obsDepth(spliceGraph, graph, exonDepthFile, junctionCountFile);

    //initialize alpha
    obsDepth.initAlpha();

    //filter cells with low expression, low quality
    map<size_t, vector<size_t>> neighborsFiltered = obsDepth.filterCells(minDepth, lowQualFrac, maxLowQual, neighbors, filteredFile);
    obsDepth.writeNumNB(numNBFile, neighborsFiltered);

    //check if no cells remaining after filtering for low expression or low quality information
    if (neighborsFiltered.empty()){
        //delete graphs
        spliceGraph->deleteGraph();
        delete(spliceGraph);
        graph->deleteGraph();
        delete(graph);

        //update progressbar
        progressbar.update(1);

        return;
    }

    //PSI solver
    PSIsolver solver(graph, &obsDepth, neighborsFiltered, wObs, wSrcSink, wFlow, wClust);

    if (VERBOSE){
        obsDepth.writeAlpha(alphaLogDir + "/init");
    }

    int iter = 0;

    while (iter < maxIter){

        //Estep
        bool converged = solver.solveEqs();

        if (VERBOSE)
            solver.writeResults(PSILogDir + "/iter" + to_string(iter));

        //test convergence
        if (converged)
            break;

        vector<VectorXd> PSI = solver.getPSI();

        //Mstep
        obsDepth.updateAlpha(PSI);

        if (VERBOSE)
            obsDepth.writeAlpha(alphaLogDir + "/iter" + to_string(iter));

        iter++;
    }

    //log nr of iterations
    if (VERBOSE)
        nIterLog.writeLine(graph->getGeneID() + "\t" + to_string(iter));

    //write PSI results to file
    solver.writeResults(outFile);

    //write alpha from last iteration to file
    obsDepth.writeAlpha(alphaFile);

    //delete graphs
    spliceGraph->deleteGraph();
    delete(spliceGraph);
    graph->deleteGraph();
    delete(graph);

    //update progressbar
    progressbar.update(1);
}

CompareArgs::CompareArgs(Progressbar &progressbar, string trueGraphFile, string truePSIFile, string estGraphFile,
                         string estPSIFile, string outDir) : Args(progressbar),
                         trueGraphFile(move(trueGraphFile)), truePSIFile(move(truePSIFile)),
                         estGraphFile(move(estGraphFile)), estPSIFile(move(estPSIFile)), outDir(move(outDir)) {}

void CompareArgs::doWork() {

    //read graphs
    SpliceGraph trueGraph(trueGraphFile);
    SpliceGraph estGraph(estGraphFile);

    //merge true splicegraph and simplified splicegrap to obtain a merged graph that contains all nodes and edges of both graphs
    SpliceGraph mergedGraph = mergeGraphs(trueGraph, estGraph);

    //write true & computed flow to file
    ofstream trueOFS(outDir + "/truePSI");
    ofstream estOFS(outDir + "/estPSI");

    //write headers
    vector<SGNode*> orderedNodes; vector<SGEdge*> orderedEdges;
    trueOFS << "cell" << mergedGraph.getVarNames(orderedNodes, orderedEdges) << endl;
    estOFS << "cell" << mergedGraph.getVarNames(orderedNodes, orderedEdges) << endl;

    //make variable mapping
    map<string, string> trueMap = mergedGraph.getVarMap(trueGraph);
    map<string, string> estMap = mergedGraph.getVarMap(estGraph);

    //remove graphs from memory
    trueGraph.deleteGraph();
    estGraph.deleteGraph();

    //read true cell flow
    PSIReader trueReader(truePSIFile);
    map<string, map<string, double>> truePSI = trueReader.readAll();

    vector<string> trueCells;
    for (auto const &it : truePSI)
        trueCells.push_back(it.first);

    //read CRF cell flow
    PSIReader estReader(estPSIFile);
    map<string, map<string, double>> estPSI = estReader.readAll();

    vector<string> estCells;
    for (auto const &it : estPSI)
        estCells.push_back(it.first);

    vector<vector<string>> cellsOverlap = Util::getNonOverlappingSets(trueCells, estCells);

    //write skipped cells to file
    ofstream skippedOFS(outDir + "/skippedCells");
    for (auto const &cell : cellsOverlap[0])
        skippedOFS << cell << endl;
    skippedOFS.close();

    double avError = 0, avAbsError = 0;
    double avNodeError = 0, avAbsNodeError = 0;
    double avEdgeError = 0, avAbsEdgeError = 0;
    double avObsEdgeError = 0, avAbsObsEdgeError = 0;
    for (auto const &cell : cellsOverlap[2]){

        //write cellName
        trueOFS << cell;
        estOFS << cell;

        //write node PSI values
        for (auto const &node : orderedNodes){
            string nodeName = node->getName();
            double truePSIVal = truePSI[cell][trueMap[nodeName]];
            double estPSIVal = estPSI[cell][estMap[nodeName]];

            //write to file
            trueOFS << "\t" << truePSIVal;
            estOFS << "\t" << estPSIVal;

            //compute error
            double error = truePSIVal - estPSIVal;
            avError += error;
            avAbsError += abs(error);
            avNodeError += error * node->getLength();
            avAbsNodeError += abs(error) * node->getLength();
        }

        //write edge PSI values
        for (auto const &edge : orderedEdges){
            string edgeName = edge->getName();
            double truePSIVal = truePSI[cell][trueMap[edgeName]];
            double estPSIVal = estPSI[cell][estMap[edgeName]];

            //write to file
            trueOFS << "\t" << truePSIVal;
            estOFS << "\t" << estPSIVal;

            //compute error
            double error = truePSIVal - estPSIVal;
            avError += error;
            avAbsError += abs(error);
            avEdgeError += error;
            avAbsEdgeError += abs(error);

            if (not edge->isSourceSink()){
                avObsEdgeError += error;
                avAbsObsEdgeError += abs(error);
            }
        }

        trueOFS << endl;
        estOFS << endl;
    }

    trueOFS.close();
    estOFS.close();

    //get sum of all node lengths
    size_t sumNodeLen = 0;
    for (auto node : orderedNodes){
        sumNodeLen += node->getLength();
    }

    //write average and average absolute error to file
    size_t nCellsComputed = cellsOverlap[2].size();
    ofstream errorOFS(outDir + "/averageError");
    errorOFS << "avError" << "\t" << avError / double(nCellsComputed * mergedGraph.getNumVars()) << endl;
    errorOFS << "avAbsError" "\t" << avAbsError / double(nCellsComputed * mergedGraph.getNumVars()) << endl;
    errorOFS << "avNodeError" << "\t" << avNodeError / double(nCellsComputed * sumNodeLen) << endl;
    errorOFS << "avAbsNodeError" << "\t" << avAbsNodeError / double(nCellsComputed * sumNodeLen) << endl;
    errorOFS << "avEdgeError" << "\t" << avEdgeError / double(nCellsComputed * mergedGraph.getNumEdges()) << endl;
    errorOFS << "avAbsEdgeError" << "\t" << avAbsEdgeError / double(nCellsComputed * mergedGraph.getNumEdges()) << endl;
    errorOFS << "avObsEdgeError" << '\t' << avObsEdgeError / double(nCellsComputed * mergedGraph.getNumObsEdges()) << endl;
    errorOFS << "avAbsObsEdgeError" << '\t' << avAbsObsEdgeError / double(nCellsComputed * mergedGraph.getNumObsEdges()) << endl;
    errorOFS.close();

    //cleanup merged graph
    mergedGraph.deleteGraph();

    progressbar.update(1);

}

//------------------------------//
// implementation of threadpool //
//------------------------------//

Pool::Pool() {
    //lock mutex
    unique_lock<mutex> lock(pool_mutex);

    //set value finished
    finished = false;
}

void Pool::addJob(Args* arg) {
    //lock mutex
    unique_lock<mutex> lock(pool_mutex);

    //add job to queue
    jobQueue.push(arg);

    //remove lock
    lock.unlock();

    //notify one waiting thread that new job is available
    condition.notify_one();
}

void Pool::doWork() {
    while(true){
        //lock mutex
        unique_lock<mutex> lock(pool_mutex);

        //wait for available jobs or until pool is finished
        condition.wait(lock, [this]() {return !jobQueue.empty() || finished;});

        if(finished && jobQueue.empty()){
            return; //finish thread, lock released automatically
        }

        //get job from queue
        Args* arg = jobQueue.front();
        jobQueue.pop();

        //remove lock
        lock.unlock();

        //do work
        arg->doWork();
    }
}

void Pool::terminatePool() {
    //lock mutex
    unique_lock<mutex> lock(pool_mutex);

    //set value finished
    finished = true;

    //remove lock
    lock.unlock();

    //notify all waiting threads
    condition.notify_all();
}


//-----------------//
// Utility classes //
//-----------------//

SJReader::SJReader(const string& SJFile) {
    //open SJ count file
    ifs = fstream(SJFile.c_str());
    if (!ifs)
    throw ios_base::failure("Cannot open file " + SJFile);
}

void SJReader::close(){
    ifs.close();
}

bool SJReader::readLine() {
    string line;
    string tmp;
    if (getline(ifs, line)){
        //read line
        istringstream iss(line);
        iss >> chrName >> start >> end >> strand >> tmp >> annot >> count;

        //convert strand
        if (strand == "1"){ //FWD strand
            strand = "+";
        }else if (strand == "2"){ //REV strand
            strand = "-";
        }else{ //UNK strand
            strand = ".";
        }

        junctionName = chrName + "_" + to_string(start) + "_" + to_string(end);
        return true;
    }
    return false;
}


SJCountReader::SJCountReader(const string& SJCountFile) {
    //open SJ count file
    ifs = fstream(SJCountFile.c_str());
    if (!ifs)
        throw ios_base::failure("Cannot open file " + SJCountFile);

    //skip first line containing cellIDs
    string tmp;
    getline(ifs, tmp);
}

void SJCountReader::close() {
    ifs.close();
}

bool SJCountReader::readLine() {
    string line;
    int count;
    if (getline(ifs, line)){
        //read line
        istringstream iss(line);
        iss >> junctionName;
        while (iss >> count)
            counts.push_back(count);
        return true;
    }
    return false;
}

PSIReader::PSIReader(const string& fileName) {
    //open file
    ifs = fstream(fileName.c_str());
    if (!ifs)
        throw ios_base::failure("Cannot open file " + fileName);

    //read header
    string line, cell;
    getline(ifs, line);

    //read tab-separated line
    vector<string> splitLine = Util::splitString(line, '\t');

#ifdef DEBUG
    //read line
    cell = splitLine[0];
    assert (cell == "cell");
#endif

    //read rest of line
    for (int i = 1; i < splitLine.size(); i++)
        varNames.push_back(splitLine[i]);
}

void PSIReader::close() {
    ifs.close();
}

bool PSIReader::readLine() {

    string line, cellName;

    if (getline(ifs, line)){

        //read tab-separated line
        vector<string> splitLine = Util::splitString(line, '\t');

        //read line
        cellName = splitLine[0];

        //add 0 PSI for variable that does not exist (varName == "")
        PSI[cellName][""] = 0;

        //add PSI value for each node/edge
        for (int i = 1; i < splitLine.size(); i++){
            PSI[cellName][varNames[i-1]] = stod(splitLine[i]);
        }

        return true;
    }
    return false;
}

map<string, map<string, double>> PSIReader::readAll() {

    //read line by line
    while (readLine())
        continue;

    close();

    return PSI;
}