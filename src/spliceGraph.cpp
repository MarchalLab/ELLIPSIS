
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <queue>

#include "spliceGraph.h"
#include "util.h"

using namespace std;

SpliceGraph::SpliceGraph(string chrName, size_t geneStart, size_t geneEnd, string geneID, string geneName, size_t strandInt,
                         size_t longestTranscriptLen, size_t shortestTranscriptLen) :
        chrName(move(chrName)), geneID(move(geneID)), geneName(move(geneName)), geneStart(geneStart), geneEnd(geneEnd),
        longestTranscriptLen(longestTranscriptLen), shortestTranscriptLen(shortestTranscriptLen) {

    //Strand can only have values: 0(undefined), 1(+), 2(-)
    assert(strandInt >= 0 && strandInt <= 2);

    //convert integer strand to character strand
    string strandString;
    if (strandInt == 0){
        strandString = '.';
    } else if (strandInt == 1){
        strandString = '+';
    } else{
        strandString = "-";
    }

    //genestart should be smaller than geneEnd
    assert(geneStart < geneEnd);

    sourceNode = new SourceSGNode(strand);
    sinkNode = new SinkSGNode(strand);
    addNode(sourceNode);
    addNode(sinkNode);

}

SpliceGraph::SpliceGraph(string chrName, size_t geneStart, size_t geneEnd, string geneID, string geneName, const string& strand,
                         size_t longestTranscriptLen, size_t shortestTranscriptLen) :
        chrName(move(chrName)), geneID(move(geneID)), geneName(move(geneName)), geneStart(geneStart), geneEnd(geneEnd), strand(strand),
        longestTranscriptLen(longestTranscriptLen), shortestTranscriptLen(shortestTranscriptLen) {

    //genestart should be smaller than geneEnd
    assert(geneStart < geneEnd);

    assert (strand == "." | strand == "+"  | strand== "-" );

    sourceNode = new SourceSGNode(strand);
    sinkNode = new SinkSGNode(strand);
    addNode(sourceNode);
    addNode(sinkNode);
}

SpliceGraph::SpliceGraph(const string& file){
    fstream ifs(file.c_str());
    if (!ifs)
        throw ios_base::failure("Cannot open file " + file);
    string line;

    //first line contains general info
    getline(ifs, line);
    istringstream issHeader(line);
    issHeader >> chrName >> geneID >> geneName >> strand >> geneStart >> geneEnd >> longestTranscriptLen >> shortestTranscriptLen;

    //add source and sink nodes
    sourceNode = new SourceSGNode(strand);
    sinkNode = new SinkSGNode(strand);
    addNode(sourceNode);
    addNode(sinkNode);

    //other lines contain edges
    string edgeStart, edgeEnd, tmp;
    while(getline(ifs, line)){
        istringstream iss(line);
        iss >> edgeStart >> tmp >> edgeEnd;

        //get start node of edge and add if it does not yet exist
        SGNode* startNode = getNode(edgeStart);
        if (startNode == nullptr){
            startNode = new SGNode(edgeStart);
            addNode(startNode);
        }

        //get end node of edge and add if it does not yet exist
        SGNode* endNode = getNode(edgeEnd);
        if (endNode == nullptr){
            endNode = new SGNode(edgeEnd);
            addNode(endNode);
        }

        //add new edge
        SGEdge* newEdge = new SGEdge(startNode, endNode, strand);
        edges.push_back(newEdge);
    }

}

SpliceGraph::SpliceGraph() : strand("."), longestTranscriptLen(0),
                             shortestTranscriptLen(0), geneStart(0), geneEnd(0) {
    sourceNode = nullptr;
    sinkNode = nullptr;
}

SpliceGraph::SpliceGraph(const SpliceGraph& g) : SpliceGraph(g.chrName, g.geneStart, g.geneEnd,
                                                             g.geneID, g.geneName,g.strand,
                                                             g.longestTranscriptLen,
                                                             g.shortestTranscriptLen){

    //mapping of old nodes to new nodes
    map<SGNode*, SGNode*> nodeMap;

    for (auto const &node : g.nodes){
        SGNode* newNode;
        if (node->isSource()){
            newNode = this->sourceNode;
        }else if (node->isSink()){
            newNode = this->sinkNode;
        }else{
            newNode = new SGNode(*node);
            this->addNode(newNode);
        }
        nodeMap.insert({node, newNode});
    }

    for (auto const &edge : g.edges){
        SGNode* newStartNode = nodeMap.find(edge->getStartNode())->second;
        SGNode* newEndNode = nodeMap.find(edge->getEndNode())->second;

        SGEdge* newEdge = new SGEdge(newStartNode, newEndNode, strand);
        this->edges.push_back(newEdge);
    }
}

void SpliceGraph::addNode(SGNode* node){

#ifdef DEBUG
    //make sure node is not already in graph
    assert(find(nodes.begin(), nodes.end(), node) == nodes.end());
    if (! node->isSource() & ! node->isSink())
        assert(node->getStartPos() != 0 & node->getStartPos() != std::numeric_limits<size_t>::max());
#endif

    nodes.push_back(node);
}

SGNode* SpliceGraph::addEdge(SGNode* prev, SGNode* newNode){

#ifdef DEBUG
    assert(find(nodes.begin(), nodes.end(), prev) != std::end(nodes)); //previous node should exist when adding edge to splicegraph"

    if (strand == "+") assert(compareNodes(prev, newNode));  //forward strand: previous node should be smaller than new node
    if (strand == "-") assert(!compareNodes(prev, newNode)); //reverse strand: previous node should be larger than new node

#endif

    for(auto it = nodes.begin(); it != nodes.end(); it ++){
        SGNode* thisNode = *it;
        //check if there is an overlap between an existing exon and the new exon
        if (thisNode->overlaps(newNode)){
            //check if new node starts before existing node
            if (newNode->startsBefore(thisNode)){

                // split newNode into 2 distinct nodes
                SGNode* E1;
                if (strand == "+"){ //forward strand
                    E1 = new SGNode(newNode->getStartPos(), thisNode->getStartPos()-1); //1st partial node for newNode
                    newNode->setStartPos(thisNode->getStartPos()); //2nd partial node for newNode
                }else{ //reverse strand
                    E1 = new SGNode(thisNode->getStartPos(), newNode->getEndPos());
                    newNode->setEndPos(thisNode->getStartPos()-1);
                }

                //add edge to first part of newNode
                E1 = this->addEdge(prev, E1);

                //add edge from E1 to rest of newNode
                newNode = this->addEdge(E1, newNode);

            }
            //check if existing node starts before new node
            else if (thisNode->startsBefore(newNode)){

                // split the existing node into 2 distinct nodes
                vector<SGNode*> splitNodes = splitNode(thisNode, newNode->getStartPos()); //split the existing (it) node

                //remove existing node from the graph
                nodes.erase(it);
                delete(thisNode);

                //replace the existing node with the 2 partial nodes
                addNode(splitNodes[0]);
                addNode(splitNodes[1]);

                //add edge between new partial nodes
                SGEdge* newEdge = new SGEdge(splitNodes[0], splitNodes[1], strand);
                edges.push_back(newEdge);

                //add the same newNode to the modified graph
                newNode = this->addEdge(prev, newNode);

            }
            //existing node ends before new node
            else if(thisNode->endsBefore(newNode)){

                // split newNode into 2 distinct nodes and add them both to the graph
                SGNode* E1;
                if (strand == "+"){ //forward strand
                    E1 = new SGNode(newNode->getStartPos(), thisNode->getEndPos()); //1st partial node for newNode
                    newNode->setStartPos(thisNode->getEndPos() + 1); //2nd partial node for newNode
                }else{ //reverse strand
                    E1 = new SGNode(thisNode->getEndPos() + 1, newNode->getEndPos());
                    newNode->setEndPos(thisNode->getEndPos());
                }

                //add edge to first part of newNode
                E1 = this->addEdge(prev, E1);

                //add edge from E1 to rest of newNode
                newNode = this->addEdge(E1, newNode);

            }
            //new node ends before existing node
            else if(newNode->endsBefore(thisNode)){

                //split existing node into 2 distinct nodes
                vector<SGNode*> splitNodes = splitNode( thisNode, newNode->getEndPos()+1); //split the existing (it) node into 2 distinct nodes

                //remove existing node from the graph
                nodes.erase(it);
                delete(thisNode);

                //replace the existing node with the 2 partial nodes
                addNode(splitNodes[0]);
                addNode(splitNodes[1]);

                //add edge between new partial nodes
                SGEdge* newEdge = new SGEdge(splitNodes[0], splitNodes[1], strand);
                edges.push_back(newEdge);

                //add the same newNode to the modified graph
                newNode = this->addEdge(prev, newNode);

            }
            //new node is exactly equal to existing node
            else{
                if (newNode != sourceNode & newNode != sinkNode){
                    //remove new node
                    delete(newNode);
                    newNode = thisNode;
                }
                if (!edgeExists(prev, thisNode)){
                    //create new edge
                    SGEdge* newEdge = new SGEdge(prev, newNode, strand);
                    edges.push_back(newEdge);
                }

            }
            return newNode;
        }
    }

    //no overlap found
    addNode(newNode);
    SGEdge* newEdge = new SGEdge(prev, newNode, strand);
    edges.push_back(newEdge);

    return newNode;
}

SGEdge* SpliceGraph::addEdgeBetweenExistingNodes(SGNode* node1, SGNode* node2){
#ifdef DEBUG
    auto it1 = std::find(nodes.begin(), nodes.end(), node1);
    auto it2 = std::find(nodes.begin(), nodes.end(), node2);
    assert(it1 != std::end(nodes) && it2 != std::end(nodes)); //nodes should exist

    if (strand == "+") assert(compareNodes(node1, node2));  //forward strand: previous node should be smaller than new node
    if (strand == "-") assert(!compareNodes(node1, node2)); //reverse strand: previous node should be larger than new node

#endif

    SGEdge* newEdge = new SGEdge(node1, node2, strand);
    edges.push_back(newEdge);
    return newEdge;
}

std::ostream &operator<<(std::ostream &out, SpliceGraph& sg)
{
    //general info
    out << sg.chrName << "\t" << sg.geneID << "\t" << sg.geneName << "\t" << sg.strand << "\t" << sg.geneStart <<
        "\t" << sg.geneEnd << "\t" << sg.longestTranscriptLen << "\t" << sg.shortestTranscriptLen << endl;

    //edges
    sg.orderEdges();
    for (auto const &edge : sg.edges){
        out << *edge << endl;
    }

    return out;
}

bool SpliceGraph::isEmpty() {
    return nodes.empty();
}

SourceSGNode* SpliceGraph::getSourceNode() const{
    return sourceNode;
}

SinkSGNode* SpliceGraph::getSinkNode() const{
    return sinkNode;
}

string SpliceGraph::getChrName() const{
    return chrName;
}

string SpliceGraph::getStrand() const{
    return strand;
}

size_t SpliceGraph::getGeneStart() const{
    return geneStart;
}

size_t SpliceGraph::getGeneEnd() const{
    return geneEnd;
}

size_t SpliceGraph::getLongestTranscriptLen() const{
    return longestTranscriptLen;
}

void SpliceGraph::setLongestTranscriptLen(size_t val) {
    longestTranscriptLen = val;
}

size_t SpliceGraph::getShortestTranscriptLen() const{
    return shortestTranscriptLen;
}

void SpliceGraph::setShortestTranscriptLen(size_t val){
    shortestTranscriptLen = val;
}

string SpliceGraph::getGeneID() const{
    return geneID;
}

vector<SGEdge*> SpliceGraph::getEdges() const{
    return edges;
}

vector<SGNode*> SpliceGraph::getNodes() const{
    return nodes;
}

size_t SpliceGraph::getNumNodes() const{
    if (nodes.empty())
        return 0;
    return nodes.size();
}

size_t SpliceGraph::getNumEdges() const{
    if (edges.empty()){
        return 0;
    }
    return edges.size();
}

size_t SpliceGraph::getNumObsEdges() const{
    size_t numSrcEdges = sourceNode->getOutEdges().size();
    size_t numSinkEdges = sinkNode->getOutEdges().size();
    return getNumEdges() - (numSrcEdges + numSinkEdges);
}

size_t SpliceGraph::getNumVars() const{
    return getNumNodes() + getNumEdges();
}

void SpliceGraph::orderNodes() {
    std::sort(nodes.begin(), nodes.end(), compareNodes);
    if (strand == "-")
        std::reverse(nodes.begin(), nodes.end());
}

void SpliceGraph::orderEdges(){
    std::sort(edges.begin(), edges.end(), compareEdges);
}

bool SpliceGraph::edgeExists(SGNode* startNode, SGNode* endNode){
    for(auto const &edge : edges){
        if (edge->getStartNode() == startNode & edge->getEndNode() == endNode)
            return true;
    }
    return false;
}

bool SpliceGraph::nodeExists(size_t startPos, size_t endPos) const{
    return (getNode(startPos, endPos) != nullptr);
}

SGNode* SpliceGraph::getNode(size_t startPos, size_t endPos) const{
    for (auto const &node : nodes){
        if (node->getStartPos() == startPos & node->getEndPos() == endPos)
            return node;
    }
    return nullptr;
}

SGNode* SpliceGraph::getNode(const string& nodeName) const{
    if (nodeName == "source")
        return sourceNode;
    if (nodeName == "sink")
        return sinkNode;

    vector<string> splitNode = Util::splitString(nodeName, '-');
    return getNode(stoi(splitNode[0]), stoi(splitNode[1]));
}

SGNode* SpliceGraph::getNodePartial(const string& nodeName) const{
    if (nodeName == "source")
        return sourceNode;
    if (nodeName == "sink")
        return sinkNode;

    vector<string> splitNode = Util::splitString(nodeName, '-');
    size_t startPos = stoi(splitNode[0]);
    size_t endPos = stoi(splitNode[1]);

    //get node to which exon (partially) belongs
    for (auto const &node : nodes) {
        if (startPos >= node->getStartPos() & endPos <= node->getEndPos())
            return node;
    }
    return nullptr;

}

SGEdge* SpliceGraph::getEdge(size_t smallPos, size_t largePos) const{
    for (auto const &edge : edges){
        if (edge->getEdgeSmallPos() == smallPos && edge->getEdgeLargePos() == largePos)
            return edge;

    }
    return nullptr;
}

SGEdge* SpliceGraph::getEdge(const string& edgeName) const{
    vector<string> splitEdge = Util::splitString(edgeName, '_');

    if (splitEdge[0] == "source" | splitEdge[1] == "source"){
        for (auto const &edge : sourceNode->getOutEdges()){
            if (strand == "+")
                if(edge->getEdgeLargePos() == stoi(splitEdge[1]))
                    return edge;
            if (strand == "-")
                if (edge->getEdgeSmallPos() == stoi(splitEdge[0]))
                    return edge;
        }
        return nullptr;
    } else if (splitEdge[0] == "sink" | splitEdge[1] == "sink"){
        for (auto const &edge : sinkNode->getInEdges()){
            if (strand == "+")
                if (edge->getEdgeSmallPos() == stoi(splitEdge[0]))
                    return edge;
            if (strand == "-")
                if (edge->getEdgeLargePos() == stoi(splitEdge[1]))
                    return edge;
        }
        return nullptr;
    }

    return getEdge(stoi(splitEdge[0]), stoi(splitEdge[1]));
}

void SpliceGraph::deleteGraph() {

    //delete edges
    for (auto it : edges)
        delete(it);
    edges.clear();

    //delete nodes
    for (auto it : nodes)
        delete(it);
    nodes.clear();

}

vector<SGNode*> SpliceGraph::splitNode(SGNode* node, size_t splitPos) {

    //create 2 new partial nodes
    SGNode* E1; SGNode* E2;
    if (strand == "+"){ //forward strand
        E1 = new SGNode(node->getStartPos(), splitPos-1);
        E2 = new SGNode(splitPos, node->getEndPos());
    }else{ //reverse strand
        E1 = new SGNode(splitPos, node->getEndPos());
        E2 = new SGNode(node->getStartPos(), splitPos-1);
    }

    //notify edges of their new start/endNodes
    vector<SGEdge*> inEdges = node->getInEdges();
    vector<SGEdge*> outEdges = node->getOutEdges();
    for (auto const &inEdge : inEdges){
        inEdge->setEndNode(E1);
        E1->addInEdge(inEdge);
    }
    for (auto const &outEdge : outEdges){
        outEdge->setStartNode(E2);
        E2->addOutEdge(outEdge);
    }

    //return newly created nodes
    vector<SGNode*> v = {E1,E2};
    return v;
}

bool SpliceGraph::addJunction(const string& junctionName){

    //get small & large position
    vector<string> junctionStr = Util::splitString(junctionName, '_');
    size_t smallPos = stoi(junctionStr[1]);
    size_t largePos = stoi(junctionStr[2]);

    //check if edge exists
    for (auto const &edge : edges){
        //check if edge already exists
        if (edge->getEdgeSmallPos() == smallPos & edge->getEdgeLargePos() == largePos)
            return false;
    }

    //check if edge lies inside gene limits
    if (smallPos <= geneStart | largePos >= geneEnd)
        return false; //Don't add edge

    //get already existing nodes
    SGNode* smallNode = getSmallNode(smallPos);
    SGNode* largeNode = getLargeNode(largePos);

    //add edge
    if (strand == "+") { //forward
        addEdgeBetweenExistingNodes(smallNode, largeNode);
    }else{ //reverse
        addEdgeBetweenExistingNodes(largeNode, smallNode);
    }
    return true;

}

SGNode* SpliceGraph::getSmallNode(size_t pos) {

    //check if position at start of gene
    if (pos == geneStart){
        if (strand == "+"){ //forward strand
            return sourceNode;
        }else{ //reverse strand
            return sinkNode;
        }
    }

    //check if smallnode already exists
    for (auto const &node : nodes){
        if (node->getEndPos() + 1 == pos)
            return node;
    }

    //check if pos inside existing node
    for (auto it = nodes.begin(); it != nodes.end(); it++){
        SGNode* node = *it;
        if (node->getStartPos() < pos & node->getEndPos() >= pos){

            //split up node
            vector<SGNode*> splitNodes = splitNode(node, pos);

            //remove existing node from graph
            delete(node);
            nodes.erase(it);

            //add partial nodes to graph
            addNode(splitNodes[0]);
            addNode(splitNodes[1]);

            //add edge between partial nodes
            SGEdge* newEdge = new SGEdge(splitNodes[0], splitNodes[1], strand);
            edges.push_back(newEdge);

            if (strand == "+"){ //forward strand
                return splitNodes[0];
            } else{ //reverse strand
                return splitNodes[1];
            }
        }
    }

    //edge starts in between existing nodes (in intron)
    //get closest node that ends before pos
    SGNode* prevNode = sourceNode;
    if (strand == "-") prevNode = sinkNode;
    for (auto const &node : nodes){
        if (node->getEndPos() < pos & node->getEndPos() > prevNode->getEndPos()){
            prevNode = node;
        }
    }

    //create new node from end of closest node to pos
    SGNode* newNode;
    if (prevNode == sourceNode | prevNode == sinkNode){
        newNode = new SGNode(geneStart, pos-1);
    }else{
        newNode = new SGNode(prevNode->getEndPos()+1, pos-1);
    }
    addNode(newNode);

    //add edge between closest node and new (intronic) node
    SGEdge* newEdge;
    if (strand == "+"){ //forward strand
        newEdge = new SGEdge(prevNode, newNode, strand);
    }else{ //reverse strand
        newEdge = new SGEdge(newNode, prevNode, strand);
    }
    edges.push_back(newEdge);

    return newNode;
}

SGNode* SpliceGraph::getLargeNode(size_t pos){

    //check if position at end of gene
    if (pos == geneStart){
        if (strand == "+"){ //forward strand
            return sinkNode;
        }else{ //reverse strand
            return sourceNode;
        }
    }

    //check if largenode already exists
    for (auto const &node : nodes){
        if (node->getStartPos() - 1 == pos)
            return node;
    }

    //check if pos in existing node
    for (auto it = nodes.begin(); it != nodes.end(); it ++){
        SGNode* node = *it;
        if (node->getStartPos() <= pos & node->getEndPos() > pos){

            //split up node
            vector<SGNode*> splitNodes = splitNode(*it, pos+1);

            //remove existing node from graph
            delete(*it);
            nodes.erase(it);

            //add partial nodes to graph
            addNode(splitNodes[0]);
            addNode(splitNodes[1]);

            //add edge between partial nodes
            SGEdge* newEdge = new SGEdge(splitNodes[0], splitNodes[1], strand);
            edges.push_back(newEdge);

            if (strand == "+"){ //forward strand
                return splitNodes[1];
            } else{ //reverse strand
                return splitNodes[0];
            }
        }
    }

    //edge starts in between existing nodes (in intronic region)
    //get closest existing node
    SGNode* nextNode = sinkNode;
    if (strand == "-") nextNode = sourceNode;
    for (auto const &node : nodes){
        if (node->getStartPos() > pos & node->getStartPos() < nextNode->getStartPos()){
            nextNode = node;
        }
    }

    //create new node from pos to closest node
    SGNode* newNode;
    if(nextNode == sourceNode | nextNode == sinkNode){
        //create new node ending at gene end
        newNode = new SGNode(pos+1, geneEnd);
    }else{
        //create new node ending just before next node
        newNode = new SGNode(pos+1, nextNode->getStartPos()-1);
    }
    addNode(newNode);

    //add edge between new (intronic) exon and closest node
    SGEdge* newEdge;
    if (strand == "+"){ //forward strand
        newEdge = new SGEdge(newNode, nextNode, strand);
    }else { //reverse strand
        newEdge = new SGEdge(nextNode, newNode, strand);
    }
    edges.push_back(newEdge);

    return newNode;
}

void SpliceGraph::removeEdge(SGEdge* edgeToRemove){

    auto pos = find(edges.begin(), edges.end(), edgeToRemove);
    edges.erase(pos);

    (edgeToRemove->getStartNode())->removeOutEdge(edgeToRemove);
    (edgeToRemove->getEndNode())->removeInEdge(edgeToRemove);

    delete(edgeToRemove);
}

void SpliceGraph::removeNode(SGNode *nodeToRemove) {

    auto pos = find(nodes.begin(), nodes.end(), nodeToRemove);
    nodes.erase(pos);

    delete(nodeToRemove);
}

size_t SpliceGraph::mergeNodes() {

    size_t numMerged = 0;

    vector<SGEdge*> oldEdges = edges;
    for (auto const &edge : oldEdges){
        //check if edge connects 2 consecutive nodes
        if (edge->getLength() == 0){
            SGNode* node1 = edge->getStartNode();
            SGNode* node2 = edge->getEndNode();

            //check if nodes can be merged
            if (! node1->isArtificial() & ! node2->isArtificial() & node1->getOutEdges().size() == 1 & node2->getInEdges().size() == 1){

                numMerged += 1;

                //remove edge between consecutive nodes
                removeEdge(edge);

                //edit node2 such that it comprises both node1 and node2
                node2->setStartPos(min(node1->getStartPos(), node2->getStartPos()));
                node2->setEndPos(max(node1->getEndPos(), node2->getEndPos()));

                //add incoming edges of node1 to incoming edges of node2
                vector<SGEdge*> inEdges = node1->getInEdges();
                for (auto const &inEdge : inEdges){
                    node2->addInEdge(inEdge);
                    inEdge->setEndNode(node2);
                }

                //remove node1 (it is now part of node2)
                removeNode(node1);

            }
        }
    }

    return numMerged;
}

bool SpliceGraph::isSpliced() {
    for (auto const &node: nodes){
        if (!node->getOutEdges().empty() && node->getOutEdges().size() > 1)
            return true; //node can be alternatively spliced
    }
    return false; //graph cannot be alternatively spliced
}

vector<SGNode*> SpliceGraph::getNodesInInterval(size_t start, size_t end){
    vector<SGNode*> out;
    for (auto const &node : nodes){
        if (node->getStartPos() >= start && node->getEndPos() <= end){
            out.push_back(node);
        }
    }
    //order nodes in interval
    sort(out.begin(), out.end(), compareNodes);
    return out;
}

SGNode* SpliceGraph::getOverlappingNode(SGEdge* edge) const{
    assert (edge->isArtificial() & ! edge->isSourceSink());
    for (auto node : nodes){
        if (node->getStartPos() < edge->getEdgeSmallPos() and node->getEndPos() > edge->getEdgeLargePos())
            return node;
    }
    return nullptr;
}

SGEdge* SpliceGraph::getExistingEdge(SGNode* startNode, SGNode* endNode) const{
    vector<SGEdge*> outEdges = startNode->getOutEdges();
    for(auto const &outEdge : outEdges){
        if (outEdge->getEndNode() == endNode){
            return outEdge;
        }
    }
    return nullptr;
}

SpliceGraph mergeGraphs(const SpliceGraph& graph1, const SpliceGraph& graph2){

    assert(graph1.geneID == graph2.geneID);
    assert(graph1.strand == graph2.strand);

    //keep which edges are already added to mergedGraph
    set<SGEdge*> edgeAdded;

    //make copy of first graph
    SpliceGraph mergedGraph(graph1);

    queue<pair<SGNode*, SGNode*>> q;

    //add each edge of graph2 to graph1 (if it does not yet exist)
    for (auto edge : graph2.getSourceNode()->getOutEdges()){
        q.push(pair<SGNode*, SGNode*>(mergedGraph.getSourceNode(), edge->getEndNode()));
        edgeAdded.insert(edge);
    }
    while(! q.empty()) {
        pair<SGNode *, SGNode *> toAdd = q.front();
        q.pop();

        SGNode *newNode = new SGNode(toAdd.second->getStartPos(), toAdd.second->getEndPos());

        SGNode *addedNode = mergedGraph.addEdge(toAdd.first, newNode);

        for (auto edge: toAdd.second->getOutEdges()) {
            if (edgeAdded.find(edge) == edgeAdded.end()){ //only add new edge if it was not already added
                q.push(pair<SGNode *, SGNode *>(addedNode, edge->getEndNode()));
                edgeAdded.insert(edge);
            }
        }
    }

    return mergedGraph;

}

map<string, string> SpliceGraph::getVarMap(SpliceGraph& graph2) const{

    //make mapping
    map<string, string> varMap;

    // Get mapping for nodes
    for (SGNode* node1 : this->getNodes()){
        string nodeName1 = node1->getName();

        //get node in graph2 to which node1 (partially belongs)
        SGNode* node2 = graph2.getNodePartial(nodeName1);

        //nodes that don't exist in smaller graph get value ""
        if (node2 == nullptr)
            varMap[nodeName1] = "";
        else
            varMap[nodeName1] = node2->getName();

    }

    // Get mapping for edges
    for (SGEdge* edge1 : this->getEdges()){
        string edgeName1 = edge1->getName();
        string edgePos1 = edge1->getPosName();

        //check if edge exist
        SGEdge* edge2 = graph2.getEdge(edgePos1);

        if (edge2 == nullptr)
            //edge that does not exist could be an internal (consecutive) edge of a merged node in graph 2
            if (edge1->isArtificial() and !edge1->isSourceSink()){
                SGNode* overlappingNode= graph2.getOverlappingNode(edge1);
                if (overlappingNode == nullptr)
                    varMap[edgeName1] = "";
                else
                    varMap[edgeName1] = overlappingNode->getName();
            }
            else
                varMap[edgeName1] = "";
        else
            varMap[edgeName1] = edge2->getName();

    }

    return varMap;
}

void SpliceGraph::writeToFile(const string& fileName){
    ofstream OFS(fileName);
    if (!OFS)
        throw ios_base::failure("Cannot write file: " + fileName);
    OFS << *this;
    OFS.close();
}

void SpliceGraph::getObservedNodes(const string& exonDepthFile, set<SGNode*>& observedNodes) const{

    //read observed exons from exonDepth files (only exons expressed in at least 1 cell are included)
    fstream ifs(exonDepthFile);

    //if file does not exist, no node is observed
    if (!ifs) return;

    //read header
    vector<SGNode*> obsNodes;
    string header, exon;
    getline(ifs, header);
    istringstream issHeader(header);
    while (issHeader >> exon){
        obsNodes.push_back(getNode(exon));
    }

    //find observed nodes : at least 1 cell has cov 1
    string line, cell;
    double sumDepth;
    while (getline(ifs, line)){
        istringstream iss(line);
        iss >> cell;
        for (auto node : obsNodes){
            iss >> sumDepth;
            if (sumDepth / double(node->getLength()) >= 1){
                observedNodes.insert(node);
            }
        }
    }

    ifs.close();

}

void SpliceGraph::getObservedEdges(const string& junctionCountFile, set<SGEdge*>& observedEdges, set<SGNode*>& observedNodes,
                                   size_t  minJunctionCount, size_t minJunctionCells) const{

    //read observed junctions from junctionCount files (only junctions expressed in at least 1 cell are included)
    fstream ifs(junctionCountFile);

    //if file does not exist, no edge is observed
    if (!ifs) return;

    //read header
    vector<SGEdge*> allEdges;
    set<SGEdge*> ambiguousEdges;
    string header, junction;
    getline(ifs, header);
    istringstream issHeader(header);
    while (issHeader >> junction){
        vector<string> junctionPos = Util::splitString(junction, '_');
        SGEdge* edge = getEdge(stoi(junctionPos[0]), stoi(junctionPos[1]));
        allEdges.push_back(edge);

        //check if edge connects observed nodes
        if (observedNodes.find(edge->getStartNode()) == observedNodes.end() or
            observedNodes.find(edge->getEndNode()) == observedNodes.end()){
                ambiguousEdges.insert(edge);
        }else{
                observedEdges.insert(edge);
        }
    }

    if (ambiguousEdges.empty()){
        ifs.close();
        return;
    }

    //find number of cells in which junction is expressed (for ambiguous edges)
    string line, cell;
    double cov;
    map<SGEdge*, size_t> numExpressed;
    while (getline(ifs, line)){
        istringstream iss(line);
        iss >> cell;
        for (auto edge : allEdges){
            iss >> cov;
            if ((ambiguousEdges.find(edge) != ambiguousEdges.end())  & (cov >= double(minJunctionCount))){
                numExpressed[edge] += 1;
            }
        }
    }
    ifs.close();

    for (auto it : numExpressed) {
        SGEdge* edge = it.first;
        size_t numCells = it.second;

        if (numCells >= minJunctionCells){
            //edge is observed in enough cells, so also its start and end node is considered observed
            observedEdges.insert(edge);
            observedNodes.insert(edge->getStartNode());
            observedNodes.insert(edge->getEndNode());
        }
    }


}

void SpliceGraph::makeConnected(set<SGNode*>& observedNodes, set<SGEdge*>& observedEdges){

    //label artificial edges between source/sink and observed node as observed
    for (auto edge : edges){
        if ((edge->getStartNode() == sourceNode && observedNodes.find(edge->getEndNode()) != observedNodes.end()) |
            (edge->getEndNode() == sinkNode && observedNodes.find(edge->getStartNode()) != observedNodes.end())){
                observedEdges.insert(edge);
        }
    }

    //sort observed nodes
    vector<SGNode*> obsNodes(observedNodes.begin(), observedNodes.end());
    sort(obsNodes.begin(), obsNodes.end(), compareNodes);
    if(strand == "-") reverse(obsNodes.begin(), obsNodes.end());

    //make sure each observed node is connected by at least 1 incoming edge (loop through nodes from src to sink)
    for (auto node : obsNodes) {

        //skip source and sink node
        if (node == sourceNode || node == sinkNode)
            continue;

        //check if connected by observed incoming edges
        vector<SGEdge *> inEdges = node->getInEdges();
        bool inConnected = false;
        for (auto inEdge : inEdges) {
            if (observedEdges.find(inEdge) != observedEdges.end()) {
                inConnected = true;
                break;
            }
        }
        //label nodes and edges on shortest path(s) from any observed node as observed
        if (!inConnected) {
            queue<vector<SGEdge *>> paths;

            set<SGNode *> newObservedNodes;
            set<SGEdge *> newObservedEdges;

            for (auto inEdge : inEdges) {
                paths.push({inEdge});
            }
            size_t len = 0;
            while (!paths.empty()) {
                vector<SGEdge *> p = paths.front();
                paths.pop();

                //don't add nodes and edges of longer paths
                if (len > 0 && p.size() > len)
                    break;

                SGNode *startNode = p.back()->getStartNode();
                if (observedNodes.find(startNode) != observedNodes.end()) {
                    //add nodes and edges in path
                    for (auto e : p) {
                        newObservedEdges.insert(e);
                        newObservedNodes.insert(e->getStartNode());
                    }
                    len = p.size();

                } else {
                    vector<SGEdge *> startInEdges = startNode->getInEdges();
                    for (auto e: startInEdges) {
                        vector<SGEdge *> newPath = p;
                        newPath.push_back(e);
                        paths.push(newPath);
                    }

                }

            }

            observedNodes.insert(newObservedNodes.begin(), newObservedNodes.end());
            observedEdges.insert(newObservedEdges.begin(), newObservedEdges.end());

        }
    }

    //make sure each observed node is connected by at least 1 outgoing edge (loop through nodes from sink to src)
    for (auto node = obsNodes.rbegin(); node != obsNodes.rend(); node++) {

        //skip source and sink node
        if (*node == sourceNode || *node == sinkNode)
            continue;

        //check if connected by observed outgoing edges
        vector<SGEdge*> outEdges = (*node)->getOutEdges();
        bool outConnected = false;
        for (auto outEdge: outEdges){
            if (observedEdges.find(outEdge) != observedEdges.end()){
                outConnected = true;
                break;
            }
        }
        //label nodes and edges on shortest path(s) from any observed node as observed
        if (! outConnected){
            queue<vector<SGEdge*>> paths;

            set<SGNode*> newObservedNodes;
            set<SGEdge*> newObservedEdges;

            for (auto outEdge : outEdges){
                paths.push({outEdge});
            }
            size_t len = 0;
            while(! paths.empty()){
                vector<SGEdge*> p = paths.front();
                paths.pop();

                //don't add nodes and edges of longer paths
                if (len > 0 && p.size() > len)
                    break;

                SGNode* endNode = p.back()->getEndNode();
                if(observedNodes.find(endNode) != observedNodes.end()){
                    //add nodes and edges in path
                    for (auto e: p){
                        newObservedEdges.insert(e);
                        newObservedNodes.insert(e->getEndNode());
                    }
                    len = p.size();

                }else{
                    vector<SGEdge*> endOutEdges = endNode->getOutEdges();
                    for (auto e : endOutEdges){
                        vector<SGEdge*> newPath = p;
                        newPath.push_back(e);
                        paths.push(newPath);
                    }
                }
            }

            observedNodes.insert(newObservedNodes.begin(), newObservedNodes.end());
            observedEdges.insert(newObservedEdges.begin(), newObservedEdges.end());

        }

    }

}

string SpliceGraph::getVarNames(vector<SGNode*>& orderedNodes, vector<SGEdge*>& orderedEdges){

    string out;

    //write nodes
    orderNodes();
    for (auto node : nodes)
        out += "\t" + node->getName();

    //write edges
    orderEdges();
    for (auto edge : edges)
        out += "\t" + edge->getName();

    orderedNodes = nodes;
    orderedEdges = edges;
    return out;
}

uint1024_t SpliceGraph::getNumberOfPathsFromSrc(SGNode* node){

    //stop at src
    if (node->isSource()){
        numSrcPaths[node] = 1;
        return 1;
    }

    if (numSrcPaths.find(node) != numSrcPaths.end()){
        return numSrcPaths[node];
    }

    uint1024_t numPaths = 0;
    for (auto edge : node->getInEdges()){
        numPaths += getNumberOfPathsFromSrc(edge->getStartNode());
    }
    numSrcPaths[node] = numPaths;
    return numPaths;

}

VarGraph::VarGraph(SpliceGraph* spliceGraph){

    geneID = spliceGraph->getGeneID();
    longestTranscriptLen = spliceGraph->getLongestTranscriptLen();

    size_t varID = 0;
    numNonArtificialEdges = 0;

    spliceGraph->orderNodes();
    //add nodes
    for (auto node : spliceGraph->getNodes()){

        //add node to graph
        Node* newNode = new Node(varID, node->getLength(), node->isArtificial());
        nodeMap.insert({node, newNode});
        nodes.push_back(newNode);

        if (node->isSource())
            srcNode = newNode;
        if (node->isSink())
            sinkNode = newNode;

        varNames.push_back(node->getName());
        varID++;
    }

    spliceGraph->orderEdges();
    for (auto edge : spliceGraph->getEdges()){

        //get nodes corresponding with start/endNode
        Node* startNode = nodeMap[edge->getStartNode()];
        Node* endNode = nodeMap[edge->getEndNode()];

        //add edge to graph
        Edge* newEdge = new Edge(varID, startNode, endNode, edge->isSourceSink(), edge->isConsecutive());
        edgeMap.insert({edge, newEdge});
        edges.push_back(newEdge);

        varNames.push_back(edge->getName());
        varID++;

        if (! edge->isSourceSink() and ! edge->isConsecutive())
            numNonArtificialEdges++;
    }

    //find edges that span each node
    for (auto node : nodeMap){
        SGNode* sgNode = node.first;
        Node* varNode = node.second;
        if (sgNode->isSource()){
            for (auto edge : varNode->getEndVars())
                spanningEdges[varNode].push_back(edge);
        }else if (sgNode->isSink()){
            for (auto edge : varNode->getStartVars())
                spanningEdges[varNode].push_back(edge);
        }else{
            size_t nodeStart = sgNode->getStartPos();
            size_t nodeEnd = sgNode->getEndPos();
            for (auto sgEdge : spliceGraph->getEdges()){
                size_t startPos = sgEdge->getEdgeSmallPos();
                size_t endPos = sgEdge->getEdgeLargePos();

                //find edges spanning the node
                if ((startPos <= nodeStart) & (endPos >= nodeEnd))
                    spanningEdges[varNode].push_back(edgeMap[sgEdge]);
            }
        }
    }

    numVars = varID;
}

Var* VarGraph::getVar(SGEdge* edge) {
    if (edge == nullptr)
        return nullptr;
    return edgeMap[edge];
}

Var* VarGraph::getVar(const string& varName){

    for (auto var : getVars()){
        if (getVarName(var) == varName){
            return var;
        }
    }
    return nullptr;

}

vector<size_t> VarGraph::getNonArtificialVarIDs(){
    vector<size_t> nonArtificialVarIDs;
    for (Node* node : nodes){
        if (! node->isSrcSinkVar())
            nonArtificialVarIDs.push_back(node->getVarID());
    }
    for(auto const &edge : edges){
        if (! edge->isSrcSinkVar())
            nonArtificialVarIDs.push_back(edge->getVarID());
    }
    return nonArtificialVarIDs;
}

void VarGraph::deleteGraph(){

    //delete edges
    for (auto it : edges)
        delete(it);
    edges.clear();

    //delete nodes
    for (auto it : nodes)
        delete(it);
    nodes.clear();

}

vector<Var*> VarGraph::getSpanningEdges(Node* node){
    return spanningEdges[node];
}

uint1024_t VarGraph::getTotalNumberOfPaths(){
    vector<uint1024_t> numPathsFromSrc(numVars, 0);

    return getNumberOfPathsFromSrc(getSinkNode(), numPathsFromSrc);
}

uint1024_t VarGraph::getNumberOfPathsFromSrc(Node* node, vector<uint1024_t> &numPathsFromSrc){

    //check if already computed
    if (numPathsFromSrc[node->getVarID()] > 0){
        return numPathsFromSrc[node->getVarID()];
    }

    //stop at src
    if (node == srcNode){
        numPathsFromSrc[srcNode->getVarID()] = 1;
        return 1;
    }
    //compute number of paths
    uint1024_t numPaths = 0;
    for (auto edge : node->getStartVars()){
        vector<Var*> inNode = edge->getStartVars();
        assert(inNode.size() == 1);
        numPaths += getNumberOfPathsFromSrc(dynamic_cast<Node *>(inNode[0]), numPathsFromSrc);
    }

    //save result
    numPathsFromSrc[node->getVarID()] = numPaths;
    return numPaths;

}