

#include "var.h"
#include <cassert>


void Node::addInEdge(Var* edge){
    assert(!edge->isNode());
    startVars.push_back(edge);
}

void Node::addOutEdge(Var* edge){
    assert(!edge->isNode());
    endVars.push_back(edge);
}

Edge::Edge(size_t varID, Node* startNode, Node* endNode, bool isSrcSink, bool isConsecutive) : Var(varID, isSrcSink),
        isConsecutive(isConsecutive) {

    startVars.push_back(startNode);
    endVars.push_back(endNode);

    //add edge to start/stop nodes
    startNode->addOutEdge(this);
    endNode->addInEdge(this);
}