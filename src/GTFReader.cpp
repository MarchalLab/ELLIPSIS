
#include "GTFReader.h"

#include "spliceGraph.h"
#include "util.h"
#include <iostream>
#include <algorithm>

using namespace std;

void GTFReader::writeGraphs(const string& outDir, vector<string> genes, const string& geneInfoFile) {

    fstream ifs(gtfFile.c_str());
    if (!ifs)
        throw ios_base::failure("Cannot open file " + gtfFile);

    ofstream geneInfoOFS;
    if (! geneInfoFile.empty())
        geneInfoOFS = ofstream(geneInfoFile.c_str());

    string line;
    SGNode* prev = nullptr; //previous node in current transcript
    SpliceGraph* graph = nullptr; //current splice graph
    size_t longestTranscriptLen = 0;
    size_t shortestTranscriptLen = INT_MAX;
    size_t thisTranscriptLen = 0;

    while(getline(ifs, line)) {
        if (line.front() == '#') { //comment line
            continue;
        }
        if (line.find("\tgene\t") != string::npos) { //gene info
            if (graph != nullptr){
                //finish off last transcript from previous gene
                graph->addEdge(prev, graph->getSinkNode());
                longestTranscriptLen = max(longestTranscriptLen, thisTranscriptLen);
                shortestTranscriptLen = min(shortestTranscriptLen, thisTranscriptLen);

                graph->setLongestTranscriptLen(longestTranscriptLen);
                graph->setShortestTranscriptLen(shortestTranscriptLen);

                //write graph to file
                graph->writeToFile(outDir + "/" + graph->getGeneID());

                //delete graph
                graph->deleteGraph();
                delete(graph);
            }

            //get new gene info
            readGeneLine(line);

            if (! geneInfoFile.empty())
                geneInfoOFS <<  geneID << '\t' << chrName << '\t' << geneStart << '\t' << geneEnd << '\t' << strand << endl;

            //reset longest transcript length
            longestTranscriptLen = 0;
            shortestTranscriptLen = INT_MAX;

            //Only consider genes in geneList
            if (genes.empty()  | find(genes.begin(), genes.end(), geneID) != genes.end()){
                graph = new SpliceGraph(chrName, geneStart, geneEnd, geneID, geneName, strand, 0, 0);

                //set prev node to source node of new graph
                prev = graph->getSourceNode();
            }else{
                graph = nullptr;
                prev = nullptr;
            }



        } else if(line.find("\ttranscript\t") != string::npos){ //new transcript
            if (graph == nullptr)  continue;
            if (prev != graph->getSourceNode()){ //Arrived at new transcript

                // finish off previous transcript
                graph->addEdge(prev, graph->getSinkNode());

                //get maximum transcript length for this gene
                longestTranscriptLen = max(longestTranscriptLen, thisTranscriptLen);
                shortestTranscriptLen = min(shortestTranscriptLen, thisTranscriptLen);

                //Start new transcript
                prev = graph->getSourceNode();
            }
            thisTranscriptLen = 0;

        }else if(line.find("\texon\t") != string::npos){ //exon information
            if (graph == nullptr)  continue;
            //add node and incoming edge
            readExonLine(line);
            auto* node = new SGNode(exonStart, exonEnd);
            node = graph->addEdge(prev, node);
            thisTranscriptLen += (exonEnd - exonStart + 1);

            prev = node;
        }
    }
    if (graph != nullptr){
        //finish off last transcript of last graph
        graph->addEdge(prev, graph->getSinkNode());
        longestTranscriptLen = max(longestTranscriptLen, thisTranscriptLen);
        shortestTranscriptLen = min(shortestTranscriptLen, thisTranscriptLen);

        graph->setLongestTranscriptLen(longestTranscriptLen);
        graph->setShortestTranscriptLen(shortestTranscriptLen);

        //write graph to file
        graph->writeToFile(outDir + "/" + graph->getGeneID());

        //delete last graph
        graph->deleteGraph();
        delete(graph);
    }

    //stop reading gtf file
    ifs.close();

    if (! geneInfoFile.empty())
        geneInfoOFS.close();
}

void GTFReader::readGeneLine(const string& line){
    string tmp;

    istringstream iss(line);
    iss >> chrName >> tmp >> tmp >> geneStart >> geneEnd >> tmp >> strand;

    geneID = line.substr(line.find(geneIdField) + geneIdField.length() + 2) ;
    geneID = geneID.substr(0, geneID.find('\"'));

    geneName = line.substr(line.find(geneNameField) + geneNameField.length() + 2);
    geneName = geneName.substr(0, geneName.find('\"'));
}

void GTFReader::readTranscriptLine(const string& line){
    transcriptID = line.substr(line.find(transcriptIDField) + transcriptIDField.length() + 2) ;
    transcriptID = transcriptID.substr(0, transcriptID.find('\"'));
}

void GTFReader::readExonLine(const string& line){
    string tmp;

    stringstream iss(line);
    iss >> tmp >> tmp >> tmp >> exonStart >> exonEnd;
}

void GTFReader::getTranscriptMapping(const string& graphDir, map<SpliceGraph*, map<SGNode*, vector<string>>>& nodeMapping,
                          map<SpliceGraph*, map<SGEdge*, vector<string>>>& edgeMapping){

    //read gtf file
    fstream gtfIFS(gtfFile.c_str());
    if (!gtfIFS)
        throw ios_base::failure("Cannot open file " + gtfFile);

    string line;
    SGNode* prev = nullptr; //previous node in current transcript
    SpliceGraph* graph = nullptr; //current splice graph

    while(getline(gtfIFS, line)) {
        if (line.front() == '#') {//comment line
            continue;
        }
        if (line.find("\tgene\t") != string::npos) {//gene info
            if (graph != nullptr){

                //finish off last transcript
                //add sink
                nodeMapping[graph][graph->getSinkNode()].push_back(transcriptID);

                //add last junction
                SGEdge* lastEdge = graph->getExistingEdge(prev, graph->getSinkNode());
                edgeMapping[graph][lastEdge].push_back(transcriptID);

            }

            //get gene information
            readGeneLine(line);

            //get graph
            if (Util::fileExists(graphDir + "/" + geneID)){
                graph = new SpliceGraph(graphDir + "/" + geneID);
                prev = nullptr;
            }else{
                graph = nullptr;
            }


        } else if(line.find("\ttranscript\t") != string::npos){ //new transcript

            if (graph == nullptr)
                continue;

            if (prev != nullptr){
                //add old transcriptID to sink node
                nodeMapping[graph][graph->getSinkNode()].push_back(transcriptID);

                //add old transcriptID to last junction
                SGEdge* lastEdge = graph->getExistingEdge(prev, graph->getSinkNode());
                edgeMapping[graph][lastEdge].push_back(transcriptID);

            }

            //get new transcriptID
            readTranscriptLine(line);

            //add new transcript to source node
            nodeMapping[graph][graph->getSourceNode()].push_back(transcriptID);

            //keep reference to last node (here source node)
            prev = graph->getSourceNode();

        }else if(line.find("\texon\t") != string::npos){ //exon information

            if (graph == nullptr)
                continue;

            readExonLine(line);

            //find corresponding nodes
            vector<SGNode*> nodes = graph->getNodesInInterval(exonStart, exonEnd);
            //sort from source to sink
            if (graph->getStrand() == "-"){ //rev strand
                reverse(nodes.begin(), nodes.end());
            }

            for (auto node: nodes){
                //add transcript to node
                nodeMapping[graph][node].push_back(transcriptID);

                //add transcript to edge
                SGEdge* edge = graph->getExistingEdge(prev, node);
                edgeMapping[graph][edge].push_back(transcriptID);

                //update prev node
                prev = node;
            }

        }
    }


    if (graph != nullptr){
        //finish off last transcript of last gene
        //add sink
        nodeMapping[graph][graph->getSinkNode()].push_back(transcriptID);

        //add last junction
        SGEdge* lastEdge = graph->getExistingEdge(prev, graph->getSinkNode());
        edgeMapping[graph][lastEdge].push_back(transcriptID);
    }

    gtfIFS.close();

}

