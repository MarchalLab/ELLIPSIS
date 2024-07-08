
#include "observedDepth.h"
#include "NNLSSolver.h"
#include "util.h"
#include <list>
#include <cassert>

//Auxiliary structs
struct maxDepthObj {
    double depth;
    Node* node;
};

ObservedDepth::ObservedDepth(SpliceGraph* spliceGraph, VarGraph* graph, const string& exonDepthFile,
                             const string& junctionDepthFile) : graph(graph) {

    //get observed depth for exons
    getExonDepth(spliceGraph, exonDepthFile);

    //get observed depth for junctions
    getJunctionDepth(spliceGraph, junctionDepthFile);

}

void ObservedDepth::getExonDepth(SpliceGraph* spliceGraph, const string& exonDepthFile){

    //open file
    fstream exonIFS = fstream(exonDepthFile.c_str());
    if (!exonIFS)
        throw ios_base::failure("Cannot open file " + exonDepthFile);

    //get order of observed exons (because of partial exons : varIDs can contain duplicates)
    vector<Var*> vars;

    string line, exon;
    getline(exonIFS, line);
    istringstream iss(line);
    while (iss >> exon){
        vars.push_back(graph->getVar(spliceGraph->getNodePartial(exon)));
    }

    //read observed depth
    string cell; double d;

    while( getline(exonIFS, line)){
        istringstream issDepth(line);

        //get cell
        issDepth >> cell;
        cells.push_back(cell);

        //initialize depth for all exons AND JUNCTIONS to 0
        vector<double> thisDepth(graph->getNumVars(), 0.0);

        //get depth for each (partial) exon
        for (Var* var : vars){
            issDepth >> d;

            //exon is not included in simplified graph
            if (var == nullptr)
                continue;

            thisDepth[var->getVarID()] += d;
        }

        //divide depth by exon length
        for (Node* node : graph->getNodes()){
            size_t nodeID = node->getVarID();
            if (thisDepth[nodeID] > 0)
                thisDepth[nodeID] /= node->getLength();
        }

        cellVarDepth.push_back(thisDepth);

    }

    //close file
    exonIFS.close();
}

void ObservedDepth::getJunctionDepth(SpliceGraph* spliceGraph, const string& junctionDepthFile) {

    //open file
    fstream ifs(junctionDepthFile);

    //check if junctions are observed
    if (!ifs)
        return; //depth for all junctions already initialized to 0 in getExonDepth

    //read header
    vector<Var*> vars;

    string line, junction;
    getline(ifs, line);
    istringstream iss(line);
    while (iss >> junction){
        SGEdge* edge = spliceGraph->getEdge(junction);
        vars.push_back(graph->getVar(edge));
    }

    //read observed depths
    string cell; double d;
    while( getline(ifs, line)){
        istringstream issDepth(line);
        issDepth >> cell;

        //find cellIdx
        auto it = find(cells.begin(), cells.end(), cell);
        assert (it != cells.end());
        size_t cellIdx = it - cells.begin();

        //get depth for each junction
        for (Var* var : vars){
            issDepth >> d;

            //artificial edge between 2 merged nodes does not exist in simplified graph => skip
            if (var == nullptr)
                continue;

            cellVarDepth[cellIdx][var->getVarID()] = d ;
        }
    }

    //close file
    ifs.close();

}

void ObservedDepth::initAlpha(){

    //at least 10% of the longest transcript length needed to accurately measure maximum depth
    //if needed, take the average of the n exons with most depth until the required number of nts are reached
    size_t minNumNT = graph->getLongestTranscriptLen() / 10;

    for (auto const& depth : cellVarDepth){

        list<maxDepthObj> maxDepths;
        double minList = 0.0;
        size_t sumNT = 0;

        //get list of exons with highest depth
        for (Node* nodeVar : graph->getNodes()){

            size_t varID = nodeVar->getVarID();

            //skip source/sink nodes
            if (nodeVar->isSrcSinkVar())
                continue;

            //get depth
            double varDepth = depth[varID];

            //check if number of nucleotides already achieved and if new depth larger
            if ((sumNT < minNumNT) | (varDepth >= minList)) {
                //add node to maxDepths
                auto dIt = maxDepths.begin();

                //find position where to insert new node depth
                while (dIt != maxDepths.end() and dIt->depth > varDepth)
                    dIt++;

                //break ties
                if (dIt != maxDepths.end() and dIt->depth == varDepth){
                    if (dIt->node->getLength() > nodeVar->getLength())
                        dIt++;
                    else if (dIt->node->getLength() == nodeVar->getLength()){
                        if (dIt->node->getVarID() < varID)
                            dIt++;
                    }
                }

                //insert element
                maxDepthObj newDepthObj = {.depth = varDepth, .node = nodeVar};
                maxDepths.insert(dIt, newDepthObj);

                //recompute number of NT
                sumNT = 0;
                dIt = maxDepths.begin();
                while( (sumNT < minNumNT) & (dIt != maxDepths.end())){
                    sumNT += dIt->node->getLength();
                    dIt ++;
                }

                //remove remainder of list
                maxDepths.erase(dIt, maxDepths.end());

                //reset minList
                minList = (--maxDepths.end())->depth;

            }
        }

        //compute average depth of max depth exons
        double avgDepth = 0.0;
        for (auto d : maxDepths){
            avgDepth += d.depth * d.node->getLength();
        }
        avgDepth = avgDepth / sumNT;

        cellAlphaInit.push_back(avgDepth);
    }

    //set alpha values to initial
    cellAlphaExon = cellAlphaInit;
    cellAlphaJunction = cellAlphaInit;

}

double ObservedDepth::getAlphaExon(size_t cellIdx){
    return cellAlphaExon[cellIdx];
}

double ObservedDepth::getAlphaJunction(size_t cellIdx) {
    return cellAlphaJunction[cellIdx];
}

void ObservedDepth::updateAlpha(vector<VectorXd> &cellPSI){

    size_t nCells = cellPSI.size();
#ifdef DEBUG
    size_t nVars = graph->getNumVars();
#endif

    vector<Node*> nodes = graph->getNodes();
    vector<Edge*> edges = graph->getEdges();

    //get average PSI sums over all cross sections
    //vector<double> avCrossSection = getAvCrossSection(cellPSI);

    for (size_t cellIdx = 0; cellIdx < nCells; cellIdx++){

        //---------- get alphaExon ----------//
        //initialize AtA & Atb
        double AtA = 0;
        double Atb = 0;

        for (auto node : nodes){
            if (node->isSrcSinkVar())
                continue;

            size_t varID = node->getVarID();
            assert(varID >= 0 and varID < nVars);

            //skip extremely small exons
            if (node->getLength() <= 3)
                continue;

            //compute weights for WLS
            double weight = sqrt(node->getLength()); //use sqrt(length) of exon as weight -> more trust in larger exons

            Atb += weight * cellPSI[cellIdx](int(varID)) * weight * getExonDepth(cellIdx, varID);
            AtA += weight * cellPSI[cellIdx](int(varID)) * weight * cellPSI[cellIdx](int(varID));
        }

        cellAlphaExon[cellIdx]= (Atb/AtA);

        assert(cellAlphaExon[cellIdx] > 0);

        //---------- get alphaJunction ----------//
        //initialize AtA & Atb
        AtA = 0;
        Atb = 0;

        for (auto edge : edges){
            if (edge->isSrcSinkVar() | edge->isConsecutiveEdge())
                continue;

            size_t varID = edge->getVarID();
            assert(varID >= 0 and varID < nVars);

            //skip junctions between extremely small exons
            if (edge->getStartVars()[0]->getLength() <= 3 | edge->getEndVars()[0]->getLength() <= 3)
                continue;

            AtA += cellPSI[cellIdx](int(varID)) *  cellPSI[cellIdx](int(varID));
            Atb += cellPSI[cellIdx](int(varID)) *  getJunctionDepth(cellIdx, varID);

        }

        if (AtA == 0 | Atb == 0){

            //no junction-counts or no junctions with PSI > 0 -> use exonAlpha
            cellAlphaJunction[cellIdx] = cellAlphaExon[cellIdx];

        }else {

            cellAlphaJunction[cellIdx] = (Atb/AtA);

        }

    }
}

double ObservedDepth::getExonDepth(size_t cellIdx, size_t varID) {
    double covScaled = cellVarDepth[cellIdx][varID];
    return covScaled;
}

double ObservedDepth::getJunctionDepth(size_t cellIdx, size_t varID) {
    double covScaled = cellVarDepth[cellIdx][varID];
    return covScaled;
}

size_t ObservedDepth::getNumCells() {
    return cells.size();
}

string ObservedDepth::getCell(size_t cellIdx){
    return cells[cellIdx];
}

map<size_t, vector<size_t>> ObservedDepth::filterCells(double minDepth, map<string, double> lowQualFrac, double maxLowQual,
                                                 map<string, vector<string>>& neighbors, const string& outFile){

    map<size_t, vector<size_t>> neighborsFiltered;
    map<string, size_t> cell2idx;

    ofstream OFS(outFile);
    OFS << "cell\tinitAlpha\tlowQualFrac" << endl;

    //remove cells without enough read depth
    for (size_t cellIdx = 0; cellIdx < cells.size();){

        string cell = cells[cellIdx];

        if ((cellAlphaExon[cellIdx] < minDepth) | (neighbors.count(cells[cellIdx]) == 0) | lowQualFrac[cell] > maxLowQual){

            OFS << cell << "\t" << cellAlphaExon[cellIdx] << "\t" << lowQualFrac[cell] << endl;

            cellVarDepth.erase(cellVarDepth.begin() + int(cellIdx));
            cellAlphaExon.erase(cellAlphaExon.begin() + int(cellIdx));
            cellAlphaJunction.erase(cellAlphaJunction.begin() + int(cellIdx));
            cells.erase(cells.begin() + int(cellIdx));

        }else{
            cell2idx[cells[cellIdx]] = cellIdx;
            cellIdx++;
        }

    }

    OFS.close();

    //make cellIdx neighboring map, where too lowly expressed cells are removed
    for (const auto& it : neighbors){
        //check if cell itself is expressed
        if (! cell2idx.count(it.first))
            continue;

        size_t thisIdx = cell2idx[it.first];
        vector<size_t> thisNeighbors;
        for (const string& n : it.second){
            if (cell2idx.count(n)) //key exists
                thisNeighbors.push_back(cell2idx[n]);
        }
        neighborsFiltered[thisIdx] = thisNeighbors;
    }

    return neighborsFiltered;

}

void ObservedDepth::writeAlpha(const string& outFile){

    ofstream OFS(outFile);

    OFS << "cell\texonAlpha\tjunctionAlpha" << endl;
    for (size_t cellIdx = 0; cellIdx < cells.size(); cellIdx++)
        OFS << cells[cellIdx] << "\t" << getAlphaExon(cellIdx) << "\t" << getAlphaJunction(cellIdx) << endl;

    OFS.close();

}



void ObservedDepth::writeNumNB(const string& outFile, const map<size_t, vector<size_t>>& neighbors){
    ofstream out(outFile);

    for (const auto& it : neighbors){
        out << cells[it.first] << "\t" << it.second.size() << endl;
    }

    out.close();
}