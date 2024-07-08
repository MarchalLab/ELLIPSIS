#include <iostream>
#include <list>
#include <string>
#include <thread>
#include <map>
#include <cmath>

#include "settings.h"
#include "GTFReader.h"
#include "geneInfo.h"
#include "util.h"
#include "pool.h"
#include "SJCounter.h"


using namespace std;

//global variables
extern bool VERBOSE;

void buildRefGraphs(Settings& settings){

    cout << "\nBuilding splice graphs from reference genome\n";
    cout << "============================================\n" << endl;

    if(Util::dirExists(settings.getAnnotGraphsDir())){
        cout << "\nSplice graphs from reference genome already built.\n";
        return;
    }

    Util::startChrono();

    Util::createDir(settings.getAnnotGraphsDir());

    GTFReader gtfReader(settings.getGTFFile());

    gtfReader.writeGraphs(settings.getAnnotGraphsDir(), vector<string>() , settings.getGeneInfoFile());

    cout << "Done building from reference genome in " << Util::stopChronoStr() << endl << endl;

}

void mergeSJFiles(Settings& settings){

    cout << "\nMerging SJ files\n";
    cout << "================\n" << endl;

    if(Util::dirExists(settings.getNovelSJDir())){
        cout << "\nSJ files already merged\n";
        return;
    }

    Util::createDir(settings.getNovelSJDir());
    Util::startChrono();

    //read geneInfo
    vector<GeneInfo> genes;
    map<string, SJCounter*> countMap;
    fstream geneInfoIFS(settings.getGeneInfoFile());
    string line;
    while (getline(geneInfoIFS, line)){
        GeneInfo gene(line);
        genes.push_back(gene);
        countMap[gene.getGeneID()] = new SJCounter();
        Util::createDir(settings.getNovelSJDir() + "/" + gene.getGeneID());
    }

    //get all cells
    string SJsuffix = "SJ.out.tab";
    regex SJPattern("(.*)(SJ.out.tab)");
    vector<string> SJFiles = Util::getFileList(settings.getAlignmentDir(), SJPattern);
    size_t numCells = SJFiles.size();
    size_t chunkSize = settings.getChunkSize();
    size_t numChunks = numCells/chunkSize + (numCells % chunkSize != 0);

    //get novelJunctions added in STAR pass2
    vector<string> novelJunctions;
    if (! settings.getNovelJunctionFile().empty()){
        fstream ifs = fstream(settings.getNovelJunctionFile().c_str());
        string chrName; int start; int end;
        while (getline(ifs, line)) {
            //read line
            istringstream iss(line);
            iss >> chrName >> start >> end;

            novelJunctions.push_back(chrName + "_" + to_string(start) + "_" + to_string(end));
        }
    }

    //progressbar
    Progressbar progressbar(numCells, "Partially merging SJ files in " + to_string(numChunks) + " chunks", cout);
    progressbar.printProgress();

    //save junctions without gene and with more than 1 gene
    ThreadSafeStringSet junctionsWithoutGene;
    ThreadSafeStringSet junctionsMultGenes;

    for (size_t chunk = 0; chunk < numChunks; chunk ++){

        size_t thisChunkSize = chunkSize;
        if ((chunk == numChunks - 1) & (numCells%chunkSize != 0))
            thisChunkSize = numCells%chunkSize;

        //create 1 threadpool per chunk
        Pool pool;
        vector<thread> threadPool;
        size_t numThreads = min(settings.getNumCPU(), thisChunkSize); //avoid idle threads
        for (int i = 0; i < numThreads; i++){
            threadPool.push_back(thread(&Pool::doWork, &pool));
        }

        vector<SJArgs*> argvector;
        for (size_t i = chunk * chunkSize; i < min(numCells, (chunk + 1) * chunkSize); i++){

            //add job to threadpool for each cell in chunk
            string file = SJFiles[i];
            regex bamPattern("(" + file.substr(0, file.length() - SJsuffix.length()) + ")([a-zA-Z\\.]*)(\\.bam)");
            vector <string> bamFiles = Util::getFileList(settings.getAlignmentDir(), bamPattern);

            if (bamFiles.empty())
                throw runtime_error(file  +" has no corresponding bamfile");
            if (bamFiles.size() > 1)
                throw runtime_error(file + " has more than 1 corresponding bamfile");

            auto* args = new SJArgs(progressbar, bamFiles[0], countMap, settings.getAlignmentDir() + "/" + file,
                                      genes, junctionsWithoutGene, junctionsMultGenes, novelJunctions);
            pool.addJob(args);
            argvector.push_back(args);
        }

        //no more jobs will be added to pool
        pool.terminatePool();

        //wait for threads to finish work
        for(auto&& t : threadPool) t.join();

        //write novel SJ counts to chunk files
        for (auto const &it : countMap){
            string gene = it.first;
            string countFile = settings.getNovelSJDir() + "/" + gene + "/chunk" + to_string(chunk);
            bool novelJunctionsObserved = it.second->writeToFile(countFile);

            //reset countmap for next chunk of cells
            it.second->clear();
        }

        //delete arg objects
        for(auto it : argvector)
            delete(it);
    }

    //make sure last progress is printed
    progressbar.printProgress();
    cout << endl;

    //remove counter objects
    for (const auto& it : countMap){
        delete(it.second);
    }

    //write junctions without gene to file
    ofstream OFS_noGene(settings.getOutDir() + "/junctionsWithoutGene.tsv");
    for (auto const &junction : junctionsWithoutGene.getSet()){
        OFS_noGene << junction << endl;
    }
    OFS_noGene.close();

    //write junctions mapping to multiple genes to file
    ofstream OFS_mult(settings.getOutDir() + "/junctionsMultipleGenes.tsv");
    for (auto const &junction : junctionsMultGenes.getSet()){
        OFS_mult << junction << endl;
    }
    OFS_mult.close();

    //merge partial files
    size_t totGenes = Util::getFileList(settings.getNovelSJDir()).size();

    if (totGenes > 0) {

        //progressbar
        Progressbar progressbar2(totGenes, "Merge partial files", cout);
        progressbar2.printProgress();

        //create threadpool
        Pool pool;
        vector<thread> threadPool;
        size_t numThreads = min(settings.getNumCPU(), totGenes); //avoid idle threads
        for (int i = 0; i < numThreads; i++) {
            threadPool.push_back(thread(&Pool::doWork, &pool));
        }

        //merge partial files
        vector<SJMergeArgs *> argvector;
        Counter numGenes;
        Counter numJunctions;
        for (const auto &geneDir: Util::getFileList(settings.getNovelSJDir())) {
            auto *args = new SJMergeArgs(progressbar2, settings.getNovelSJDir() + "/" + geneDir, numGenes,
                                         numJunctions);
            pool.addJob(args);
            argvector.push_back(args);
        }

        //finish jobs
        pool.terminatePool();
        for (auto &&t: threadPool) t.join();
        for (auto it: argvector)
            delete (it);
        progressbar2.printProgress();
        cout << endl;

        cout << numJunctions.getValue() << " de novo junctions found in " << numGenes.getValue() << " genes." << endl;
    }

    cout << "\t" << junctionsMultGenes.getSize() << " novel junctions mapped to multiple genes." << endl;
    cout << "\t" << "ignored " << junctionsWithoutGene.getSize() << " observed junctions that did not map to any gene." << endl;
    cout << "Done merging SJ files in " << Util::stopChronoStr() << endl << endl;

}

void addDeNovo(Settings& settings){

    cout << "\nAdding de novo splice events\n";
    cout << "============================\n" << endl;

    if(Util::dirExists(settings.getDeNovoDir())){
        cout << "\nDe Novo edges already added\n";
        return;
    }

    Util::startChrono();

    Util::createDir(settings.getDeNovoDir());

    //get number of cells
    vector<string> cells = Util::getFileList(settings.getAlignmentDir(), regex("(.*)(SJ.out.tab)"));
    size_t nCells = cells.size();

    //get all annotated genes
    vector<string> geneNames = Util::getFileList(settings.getAnnotGraphsDir());

    size_t numGenes = geneNames.size();
    size_t numThreads = min(settings.getNumCPU(), numGenes); //avoid idle threads

    //progressbar
    Progressbar progressbar(numGenes, "Adding de novo edges", cout);

    //create threadpool
    Pool pool;
    vector<thread> threadPool;
    for (int i = 0; i < numThreads; i++){
        threadPool.push_back(thread(&Pool::doWork, &pool));
    }

    vector<DeNovoArgs*> argvector;

    progressbar.printProgress();

    Counter numDeNovo;
    Counter numGenesChanged;

    for (auto const &gene : geneNames){
        auto* args = new DeNovoArgs(progressbar, settings.getAnnotGraphsDir() + "/" + gene,
                                          settings.getNovelSJDir() + "/" + gene,
                                          settings.getDeNovoDir() + "/" + gene,
                                          settings.getMinJunctionCount(),
                                          settings.getMinJunctionCells(),
                                          numDeNovo,  numGenesChanged);
        pool.addJob(args);
        argvector.push_back(args);
    }

    //no more jobs will be added to pool
    pool.terminatePool();

    //wait for threads to finish work
    for(auto&& t : threadPool) t.join();

    //delete arg objects
    for(auto it : argvector)
        delete(it);

    progressbar.printProgress(); //make sure last progress is printed
    cout << endl;

    cout << "Added " << numDeNovo.getValue() << " new junctions to " << numGenesChanged.getValue() << " genes." << endl;
    cout << "Done adding De Novo edges in " << Util::stopChronoStr() << endl << endl;

}

void getDepth(Settings& settings){
    cout << "\nGetting read-depth\n";
    cout << "==================\n" << endl;

    if(Util::dirExists(settings.getExonDepthDir())){
        cout << "\nRead depths already generated\n";
        return;
    }

    Util::startChrono();

    regex bamPattern(".*\\.bam$");
    size_t numCells = Util::getFileList(settings.getAlignmentDir(), bamPattern).size();

    //assign reads to genes and compute depth
    string command = "python " + string(BUILDDIR) + "/getGeneDepth.py -n " + to_string(settings.getNumCPU()) +
                     " --deNovo_dir " + settings.getDeNovoDir() +
                     " --bam_dir " + settings.getAlignmentDir() +
                     " --exonDepth_dir " + settings.getExonDepthDir() +
                     " --junctionDepth_dir " + settings.getJunctionCountDir() +
                     " --geneCount_dir " + settings.getGeneCountDir() +
                     " --lowQual_dir " + settings.getLowQualDir() +
                     " --rawCountFile " + settings.getCountFile();
    if (! settings.getSelectCellsFile().empty())
        command += " --selectCellsFile " + settings.getSelectCellsFile();
    if (! settings.isPairedEnd())
        command += " --singleEnd";

    std::cout << command << std::endl;

    int ret = system(command.c_str());

    if (ret != 0)
        throw runtime_error("Depth could not be computed");

    cout << "Done getting depth in " << Util::stopChronoStr() << endl << endl;
}

void simplifyGraphs(Settings& settings){

    cout << "\nSimplifying splice graphs\n";
    cout << "=========================\n" << endl;

    string simplifiedDir = settings.getSimplifiedDir();
    if(Util::dirExists(simplifiedDir)){
        cout << "\nGraphs already simplified\n";
        return;
    }

    Util::startChrono();

    Util::createDir(simplifiedDir);

    //get geneNames
    vector<string> selectedGenes = settings.getSelectedGenes();
    vector<string> geneNames = Util::getFileList(settings.getGeneCountDir());
    if (! selectedGenes.empty()){
        for (auto it = geneNames.begin(); it != geneNames.end();){
            //only select geneNames that are in selectedGenes
            if ( std::find(selectedGenes.begin(), selectedGenes.end(), *it) == selectedGenes.end())
                it = geneNames.erase(it);
            else
                it++;
        }
    }

    size_t numGenes = geneNames.size();
    size_t numThreads = min(settings.getNumCPU(), numGenes); //avoid idle threads

    //counters
    Counter numNodesRemoved;        //number of nodes that were removed
    Counter numGraphsNodesRemoved;  //number of graphs where at least 1 node was removed
    Counter numEdgesRemoved;        //number of edges that were removed
    Counter numGraphsEdgesRemoved;  //number of graphs where at least 1 edge was removed
    Counter numMerged;              //number of merge operations
    Counter numGraphsMerged;        //number of graphs where at least 1 merge operation was performed
    Counter numNodesRemaining;      //number of remaining nodes
    Counter numEdgesRemaining;      //number of remaining edges

    vector<SimplifyArgs*> argvector;

    //progressbar
    Progressbar progressbar(numGenes, "Simplifying graphs", cout);
    progressbar.printProgress();

    //create threadpool
    Pool pool;
    vector<thread> threadPool;
    for (int i = 0; i < numThreads; i++){
        threadPool.push_back(thread(&Pool::doWork, &pool));
    }

    for (auto const& gene : geneNames){
        int i = 0;
        auto* args = new SimplifyArgs(progressbar, settings.getDeNovoDir() + "/" + gene,
                                              settings.getExonDepthDir() + "/" + gene,
                                              settings.getJunctionCountDir() + "/" + gene,
                                              settings.getGeneCountDir() + "/" + gene,
                                              settings.getMinDepth(),
                                              settings.getMinJunctionCount(),
                                              settings.getMinJunctionCells(),
                                              settings.getSimplifiedDir() + "/" + gene,
                                              numNodesRemoved, numGraphsNodesRemoved, numEdgesRemoved,
                                              numGraphsEdgesRemoved, numMerged, numGraphsMerged,
                                              numNodesRemaining, numEdgesRemaining);
        pool.addJob(args);
        argvector.push_back(args);
    }

    //no more jobs will be added to pool
    pool.terminatePool();

    //wait for threads to finish work
    for(auto&& t : threadPool) t.join();

    //delete arg objects
    for(auto it : argvector)
        delete(it);

    //make sure last progress is printed
    progressbar.printProgress();
    cout << endl;

    cout << "Removed " << numNodesRemoved.getValue() << " nodes in " << numGraphsNodesRemoved.getValue() << " genes." << endl;
    cout << "Removed " << numEdgesRemaining.getValue() << " edges in " << numGraphsEdgesRemoved.getValue() << " genes." << endl;
    cout << "Performed " << numMerged.getValue() << " merge operations in " << numGraphsMerged.getValue() << " genes." << endl;

    //get nr of graphs after simplify step
    int numGraphsRemaining = Util::getNumFiles(settings.getSimplifiedDir());

    cout << "Average number of nodes : " << numNodesRemaining.getValue() * 1.0 / numGraphsRemaining << endl;
    cout << "Average number of edges : " << numEdgesRemaining.getValue() * 1.0 / numGraphsRemaining << endl;

    cout << endl << "Done simplifying graphs in " << Util::stopChronoStr() << endl << endl;

}

void computePSI(Settings& settings){
    cout << "\nComputing PSI values\n";
    cout << "=================\n" << endl;

    if(Util::dirExists(settings.getPSIDir())){
        cout << "\nPSI values already computed\n";
        return;
    }

    Util::createDir(settings.getPSIDir());
    Util::createDir(settings.getFilteredDir());
    Util::createDir(settings.getNumNBDir());
    Util::createDir(settings.getAlphaDir());

    if (VERBOSE) {
        Util::createDir(settings.getAlphaLogDir());
        Util::createDir(settings.getPSILogDir());
    }
    Util::startChrono();

    //get neighbors
    map<string, vector<string>> neighbors = Util::readNeighborFile(settings.getNeighborFile(), settings.getAlignmentDir());

    //get total number of genes
    vector<string> geneNames = Util::getFileList(settings.getSimplifiedDir());
    vector<string> selectedGenes = settings.getSelectedGenes();
    if (! selectedGenes.empty()){
        for (auto it = geneNames.begin(); it != geneNames.end();){
            //only select geneNames that are in selectedGenes
            if ( std::find(selectedGenes.begin(), selectedGenes.end(), *it) == selectedGenes.end())
                it = geneNames.erase(it);
            else
                it++;
        }
    }
    size_t numGenes = geneNames.size();

    size_t numThreads = min(settings.getNumCPU(), numGenes); //avoid idle threads

    //define empty logwriters for non-verbose setting
    LogWriter nIterLog;   // log nr of iterations per gene
    LogWriter complexLog; //log genes ignored due to too many paths
    //logfiles for verbose setting
    if(VERBOSE) {
        nIterLog.init(settings.getIterLog());
        complexLog.init(settings.getComplexLog());
    }

    Progressbar progressbar(numGenes, "Computing PSI values", cout);
    progressbar.printProgress();

    //create threadpool
    Pool pool;
    vector<thread> threadPool;
    for (int i = 0; i < numThreads; i++){
        threadPool.push_back(thread(&Pool::doWork, &pool));
    }

    vector<PSIArgs*> argvector;
    for (auto const &gene : geneNames){

        //add job to jobpool
        auto* args = new PSIArgs(progressbar, nIterLog, complexLog,
                                    settings.getSimplifiedDir() + "/" + gene,
                                    settings.getExonDepthDir() + "/" + gene,
                                    settings.getJunctionCountDir() + "/" + gene,
                                    settings.getPSIDir() + "/" + gene,
                                    settings.getAlphaDir() + "/" + gene,
                                    settings.getAlphaLogDir() + "/" + gene,
                                    settings.getPSILogDir() + "/" + gene,
                                    settings.getFilteredDir() + "/" + gene,
                                    settings.getLowQualDir()  +"/" + gene,
                                    settings.getNumNBDir() + "/" + gene,
                                    neighbors,
                                    settings.getWObs(),
                                    settings.getWSrcSink(),
                                    settings.getWFlow(),
                                    settings.getWClust(),
                                    settings.getMinDepth(),
                                    settings.getMaxIter(),
                                    settings.getMaxLowQual(),
                                    settings.getMaxPaths());
        pool.addJob(args);
        argvector.push_back(args);

    }

    //no more jobs will be added to pool
    pool.terminatePool();

    //wait for threads to finish work
    for(auto&& t : threadPool) t.join();

    if (VERBOSE) {
        //close log writer
        nIterLog.close();
        complexLog.close();
    }

    //delete arg objects
    for(auto it: argvector)
        delete(it);

    progressbar.printProgress(); //make sure last progress is printed
    cout << endl;

    cout << "Done computing PSI values in " << Util::stopChronoStr() << endl << endl;
}

void getGroundTruth(Settings& settings){
    cout << "\n Getting ground truth\n";
    cout << "============================================\n" << endl;

    string trueGraphDir = settings.getTrueGraphsDir();
    if(Util::dirExists(trueGraphDir)){
        cout << "\nSplice graphs from true reference genome already built.\n";
        return;
    }

    Util::startChrono();

    //get true graphs from gtf file
    Util::createDir(settings.getTrueGraphsDir());
    GTFReader gtfReader(settings.getTrueGTF());
    const vector<string>& genes = settings.getSelectedGenes();
    gtfReader.writeGraphs(settings.getTrueGraphsDir(), genes, settings.getTrueGeneInfoFile());


    //get PSI values from true clusterflow
    Util::createDir(settings.getTrueCellPSI());
    Util::getTruePSI(settings.getTrueGTF(), settings.getTrueGraphsDir(),
                     settings.getTrueClusterFlow(),
                     settings.getClusterFile(),
                     settings.getTrueCellPSI());

    map<SpliceGraph*, map<SGNode*, vector<string>>> nodeMapping;
    map<SpliceGraph*, map<SGEdge*, vector<string>>> edgeMapping;
    gtfReader.getTranscriptMapping(settings.getTrueGraphsDir(), nodeMapping, edgeMapping);

    cout << "Done creating ground truth in " << Util::stopChronoStr() << endl << endl;
}

void compareGraphs(Settings& settings) { //compare true percentages with computed percentages
    cout << "\nComparing graphs \n";
    cout << "=================\n" << endl;

    if(Util::dirExists(settings.getComparisonDir())){
        cout << "\nGraphs already compared\n";
        return;
    }

    Util::createDir(settings.getComparisonDir());
    Util::startChrono();

    //get genes for which PSI is computed and true graph exists
    vector<string> geneNamesEst = Util::getFileList(settings.getPSIDir());                  //genes for which CRF is computed
    vector<string> geneNamesTrue = Util::getFileList(settings.getOutDir() + "/trueGraphs"); //genes for which true flow exists
    vector<vector<string>> nonOverlappingSets = Util::getNonOverlappingSets(geneNamesEst, geneNamesTrue);

    //write genes for which no PSI was computed
    ofstream notComputedOFS(settings.getComparisonDir() + "/genesNotComputed");
    for (const string& gene : nonOverlappingSets[1])
        notComputedOFS << gene << endl;
    notComputedOFS.close();
    cout << "Did not compute " << nonOverlappingSets[1].size() << " genes, that were in ground truth." << endl;

    //write genes for which ground truth did not exist, but PSI values were estimated
    ofstream notInTrueOFS(settings.getComparisonDir() + "/genesNotInTrue");
    for (const string& gene : nonOverlappingSets[0])
        notInTrueOFS << gene << endl;
    notInTrueOFS.close();
    cout << "Computed " << nonOverlappingSets[0].size() << " genes that were not in ground truth." << endl;

    //genes where PSI is computed and ground truth exists
    vector<string> geneNames = nonOverlappingSets[2];
    size_t numGenes = geneNames.size();

    size_t numThreads = min(settings.getNumCPU(), numGenes); //avoid idle threads

    Progressbar progressbar(numGenes, "Computing accuracy", cout);
    progressbar.printProgress();

    //create threadpool
    Pool pool;
    vector<thread> threadPool;
    for (int i = 0; i < numThreads; i++){
        threadPool.push_back(thread(&Pool::doWork, &pool));
    }

    vector<CompareArgs*> argvector;

    for (const auto& gene : geneNames){

        Util::createDir(settings.getComparisonDir() + "/" + gene);

        //add job to jobpool
        auto* args = new CompareArgs(progressbar,
                                            settings.getTrueGraphsDir() + "/" + gene,
                                            settings.getTrueCellPSI() + "/" + gene,
                                            settings.getSimplifiedDir() + "/" + gene,
                                            settings.getPSIDir() + "/" + gene,
                                            settings.getComparisonDir() + "/" + gene);
        pool.addJob(args);
        argvector.push_back(args);

    }

    //no more jobs will be added to pool
    pool.terminatePool();

    //wait for threads to finish work
    for(auto&& t : threadPool) t.join();

    //delete arg objects
    for(auto it : argvector)
        delete(it);

    //make sure last progress is printed
    progressbar.printProgress();
    cout << endl;

    cout << "Done comparing PSI values in " << Util::stopChronoStr() << endl << endl;
}

int main(int argc, char** argv) {

    try {
        Settings settings(argc, argv);

        cout << "Welcome to ELLIPSIS" ;

        #ifdef DEBUG
                cout << " (debug mode)" << endl;
        #else
                cout << " (release mode)" << endl;
        #endif

        cout << " (commit " << GIT_SHA1 << ")" << endl << endl;

        for(int i = 0; i < argc; ++i)
            cout << argv[i] << " ";
        cout << endl << endl;

        //run program
        buildRefGraphs(settings);
        mergeSJFiles(settings);
        addDeNovo(settings);
        getDepth(settings);
        simplifyGraphs(settings);
        computePSI(settings);

        //test correctness if ground truth is given
        if ( (!settings.getTrueGTF().empty()) & (!settings.getTrueClusterFlow().empty())) {
            getGroundTruth(settings);
            compareGraphs(settings);
        }

    } catch (exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    cout << "Exiting... bye!" << endl;
    return EXIT_SUCCESS;
}
