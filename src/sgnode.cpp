
#include <iostream>
#include <limits>
#include <algorithm>
#include <cassert>

#include "sgnode.h"
#include "sgedge.h"
#include "util.h"

using namespace std;

SGNode::SGNode(size_t start, size_t end) {

    assert(start >= 0 & end >= 0);
    assert(start <= end); //exon end should be larger than exon start

    startPos = start;
    endPos = end;

}

SGNode::SGNode(const string& nodeName){

    assert (nodeName != "source");
    assert (nodeName != "sink");

    vector<string> splitNode = Util::splitString(nodeName, '-');
    size_t start = stoi(splitNode[0]);
    size_t end = stoi(splitNode[1]);

    assert(start >= 0 & end >= 0);
    assert(start <= end); //exon end should be larger than exon start

    startPos = start;
    endPos = end;

}

SGNode::SGNode(const SGNode& oldNode){

    startPos = oldNode.startPos;
    endPos = oldNode.endPos;

}

size_t SGNode::getStartPos() const {
    return startPos;
}

size_t SGNode::getEndPos() const {
    return endPos;
}

void SGNode::setStartPos(size_t newStartPos){

    assert(newStartPos <= endPos); //exon end should be larger than exon start

    startPos = newStartPos;

}

void SGNode::setEndPos(size_t newEndPos){

    assert(newEndPos >= startPos); //exon end should be larger than exon start

    endPos = newEndPos;

}

vector<SGEdge*> SGNode::getInEdges(){
    return inEdges;
}
vector<SGEdge*> SGNode::getOutEdges(){
    return outEdges;
}

bool SGNode::operator<(const SGNode& node2) const{
    return (this->startPos < node2.startPos);
}

bool SGNode::operator==(const SGNode& node2) const{
    return (this->startPos == node2.startPos & this->endPos == node2.endPos);
}

bool SGNode::overlaps(SGNode* node2) const{
    return node2->getStartPos() <= endPos & node2->getEndPos() >= startPos;
}

bool SGNode::startsBefore(SGNode* node2) const{
    return startPos < node2->startPos;
}

bool SGNode::endsBefore(SGNode* node2) const{
    return endPos < node2->endPos;
}

void SGNode::addInEdge(SGEdge *newEdge) {
    if (std::find(inEdges.begin(), inEdges.end(), newEdge) == inEdges.end()){
        inEdges.push_back(newEdge);
    }else{
        cerr << "Warning: adding already existing incoming edges to node." << endl;
    }
}

void SGNode::addOutEdge(SGEdge *newEdge) {
    if (std::find(outEdges.begin(), outEdges.end(), newEdge) == outEdges.end()) {
        outEdges.push_back(newEdge);
    }else{
        cerr << "Warning: adding already existing outgoing edges to node." << endl;
    }
}

void SGNode::removeInEdge(SGEdge* edgeToRemove){
    auto pos = std::find(inEdges.begin(), inEdges.end(), edgeToRemove);
    if (pos != inEdges.end()) {
        inEdges.erase(pos);
    }else{
        cerr << "Warning: Trying to remove incoming edge that doesn't exist." << endl;
    }
}

void SGNode::removeOutEdge(SGEdge* edgeToRemove){
    auto pos = std::find(outEdges.begin(), outEdges.end(), edgeToRemove);
    if (pos != outEdges.end()) {
        outEdges.erase(pos);
    }else{
        cerr << "Warning: Trying to remove outgoing edge that doesn't exist." << endl;
    }
}

size_t SGNode::getLength(){
    return endPos - startPos + 1;
}

bool SGNode::isArtificial(){
    return false;
}

bool SGNode::isSource() {
    return false;
}

bool SGNode::isSink() {
    return false;
}

string SGNode::getName() {
    return to_string(startPos) + "-" + to_string(endPos);
}

SourceSGNode::SourceSGNode(const string& strand) : SGNode(0,0){
    if(strand == "-"){ //reverse strand
        this->startPos = std::numeric_limits<size_t>::max();
        this->endPos = std::numeric_limits<size_t>::max();
    }
}

size_t SourceSGNode::getLength() {
    return 0;
}

bool SourceSGNode::isArtificial(){
    return true;
}

bool SourceSGNode::isSource() {
    return true;
}

bool SourceSGNode::isSink() {
    return false;
}

string SourceSGNode::getName() {
    return "source";
}

SinkSGNode::SinkSGNode(const string& strand) : SGNode(0,0){
    if (strand == "+"){ //forward strand
        this->startPos = std::numeric_limits<size_t>::max();
        this->endPos = std::numeric_limits<size_t>::max();
    }
}

size_t SinkSGNode::getLength(){
    return 0;
}

bool SinkSGNode::isArtificial(){
    return true;
}

bool SinkSGNode::isSource() {
    return false;
}

bool SinkSGNode::isSink() {
    return true;
}

string SinkSGNode::getName() {
    return "sink";
}

bool compareNodes(SGNode* node1, SGNode* node2){
    return node1->getStartPos() < node2->getStartPos();
}

std::ostream &operator<<(std::ostream &out, const SGNode &n){
    out << n.getStartPos() << "-" << n.getEndPos() ;
    return out;
}