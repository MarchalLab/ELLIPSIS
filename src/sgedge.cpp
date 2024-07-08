
#include <iostream>
#include <cassert>

#include "sgedge.h"
#include "sgnode.h"

using namespace std;

SGEdge::SGEdge(SGNode* start, SGNode* end, const string& strand){
    startNode = start;
    endNode = end;
    this->strand = strand;

    startNode->addOutEdge(this);
    endNode->addInEdge(this);

}

SGNode* SGEdge::getStartNode(){
    return(startNode);
}

SGNode* SGEdge::getEndNode(){
    return(endNode);
}

string SGEdge::getStrand(){
    return(strand);
}

void SGEdge::setStartNode(SGNode* newStartNode){
    if (strand == "+"){ //forward strand
        assert(startNode->getEndPos() == newStartNode->getEndPos()); //start position of edge needs to be conserved when changing start node
    }else{ //reverse strand
        assert(startNode->getStartPos() == newStartNode->getStartPos()); //start position of edge needs to be conserved when changing start node
    }
    startNode = newStartNode;
}

void SGEdge::setEndNode(SGNode *newEndNode) {
    if (strand == "+"){ //forward strand
        assert(endNode->getStartPos() == newEndNode->getStartPos()); //end position of edge needs to be conserved when changing end node
    }else { //reverse strand
        assert(endNode->getEndPos() == newEndNode->getEndPos()); //end position of edge needs to be conserved when changing end node
    }
    endNode = newEndNode;
}

size_t SGEdge::getEdgeSmallPos() const {
    if (strand == "+"){ //forward strand
        return startNode->getEndPos() + 1;
    }else{ //reverse strand
        return endNode->getEndPos() + 1;
    }
}

size_t SGEdge::getEdgeLargePos() const {
    if (strand == "+"){ //forward strand
        return endNode->getStartPos() - 1;
    }else{
        return startNode->getStartPos() - 1;
    }
}

size_t SGEdge::getLength(){
    if (strand == "+"){ //forward strand
        return endNode->getStartPos() - startNode->getEndPos() -1;
    }else{ //reverse strand
        return startNode->getStartPos() - endNode->getEndPos() -1;
    }
}

bool SGEdge::isArtificial(){

    if (getLength() == 0)
        return true;

    if (startNode->isArtificial() || endNode->isArtificial())
        return true;

    return false;
}

bool SGEdge::isConsecutive() {
    if (getLength() == 0)
        return true;
    return false;
}

bool SGEdge::isSourceSink() const{
    if (startNode->isArtificial() || endNode->isArtificial())
        return true;

    return false;
}

string SGEdge::getName() const{
    return startNode->getName() + " -> " + endNode->getName();
}

string SGEdge::getPosName() const{
    if (startNode->isSource()){
        if  (strand == "+")
            return "source_" + to_string(getEdgeLargePos());
        else
            return to_string(getEdgeSmallPos()) + "_source";
    }
    if (endNode->isSink()){
        if (strand == "+")
            return to_string(getEdgeSmallPos()) + "_sink";
        else
            return "sink_" + to_string(getEdgeLargePos());
    }
    return to_string(getEdgeSmallPos()) + "_" + to_string(getEdgeLargePos());
}

bool SGEdge::operator==(const SGEdge& edge2){
    return (this->startNode == edge2.startNode & this->endNode == edge2.endNode);
}

bool compareEdges(SGEdge* edge1, SGEdge* edge2){
    bool ret;

    if (edge1 == nullptr | edge2 == nullptr){
        cout << "test";
    }

    if (edge1->getStartNode() == edge2->getStartNode()){ //if startnodes are the same, comparison based on end nodes
        assert(edge1->getEndNode() != edge2->getEndNode()); //if both end and start node are the same: 2x same edge in graph
        ret = compareNodes(edge1->getEndNode(), edge2->getEndNode());

    }else{
        ret = compareNodes(edge1->getStartNode(), edge2->getStartNode());
    }

    if (edge1->getStrand() == "+"){ //forward strand
        return ret;
    }else{ //reverse strand
        return !ret; //reverse strand => order from right to left
    }


}

std::ostream &operator<<(std::ostream &out, const SGEdge &e){
    out << e.getName();
    return out;
}
