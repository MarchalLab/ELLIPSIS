

#ifndef ELLIPSIS_SPLICEGRAPH_H
#define ELLIPSIS_SPLICEGRAPH_H

#include <boost/multiprecision/cpp_int.hpp>

#include <string>
#include <vector>
#include <set>

#include "sgedge.h"
#include "sgnode.h"
#include "var.h"

using namespace boost::multiprecision;
using namespace std;

class SpliceGraph {

private:
    string chrName;
    string geneID;
    string geneName;
    string strand; //+ or - strand (. for undefined)
    size_t geneStart;
    size_t geneEnd;
    size_t longestTranscriptLen;
    size_t shortestTranscriptLen;

    SourceSGNode* sourceNode;
    SinkSGNode* sinkNode;
    vector<SGNode*> nodes;
    vector<SGEdge*> edges;

    map<SGNode*, uint1024_t> numSrcPaths;

public:
    /**
     * Constructor using integer strands
     * @param geneID
     * @param geneName
     * @param strand integer values 0(undefined), 1(forward), 2(reverse)
     */
    SpliceGraph(string chrName, size_t geneStart, size_t geneEnd, string geneID, string geneName, size_t strand,
                size_t longestTranscriptLen, size_t shortestTranscriptLen);

    /**
     * Constructor using char strands
     * @param geneID
     * @param geneName
     * @param strand String values .(undefined), +(forward), -(reverse)
     */
    SpliceGraph(string chrName, size_t geneStart, size_t geneEnd, string geneID, string geneName, const string& strand,
                size_t longestTranscriptLen, size_t shortestTranscriptLen);

    /**
     * constructor using file
     * @param file
     */
    explicit SpliceGraph(const string& file);

    /**
     * Default constructor
     */
    SpliceGraph();

    /**
     * Copy constructor
     * Also makes copies of nodes and edges
     */
    SpliceGraph(const SpliceGraph& p1);

    /**
     * add node to splicegraphs
     * @param node
     */
    void addNode(SGNode* node);

    /**
     * Add newNode and edge between prev and newNode to splicegraph
     * @param prev node already in splicegraph
     * @param newNode new node to be added to splicegraph
     * @return pointer to newNode
     */
    SGNode* addEdge(SGNode* prev, SGNode* newNode);

    /**
     * Add a new edge between node1 and node2
     * @param node1
     * @param node2
     */
    SGEdge* addEdgeBetweenExistingNodes(SGNode* node1, SGNode* node2);

    /**
     * Operator<< overloading
     * @param out Output file stream
     * @param ga splicegraph to display
     * @return Output file stream
     */
    friend std::ostream &operator<<(std::ostream &out,
                                    SpliceGraph& sg);

    /**
     * Check if graph is empty
     * @return True if no nodes in graph
     */
    bool isEmpty();

    SourceSGNode* getSourceNode() const;

    SinkSGNode* getSinkNode() const;

    string getChrName() const;

    string getStrand() const;

    size_t getGeneStart() const;

    size_t getGeneEnd() const;

    size_t getLongestTranscriptLen() const;

    void setLongestTranscriptLen(size_t val);

    size_t getShortestTranscriptLen() const;

    void setShortestTranscriptLen(size_t val);

    string getGeneID() const;

    vector<SGEdge*> getEdges() const;

    vector<SGNode*> getNodes() const;

    size_t getNumNodes() const;

    size_t getNumEdges() const;

    /**
     * get number of observable edges (all edges except the ones connecting to src/sink)
     * @return
     */
    size_t getNumObsEdges() const;

    size_t getNumVars() const;

    /**
     * order nodes based on position (from source to sink)
     *      from smaller pos to larger pos for forward strand
     *      from larger pos to smaller pos for reverse strand
     */
    void orderNodes();

    /**
     * order edges: from left to right for forward strand
     *              from right to left for reverse strand
     */
    void orderEdges();

    /**
     * Check if splice graph already contains an edge between these nodes
     * @param startNode
     * @param endNode
     * @return true if the edge is already part of the graph
     */
    bool edgeExists(SGNode* startNode, SGNode* endNode);

    /**
     * Check if splice graph already contains a node with the same start and end position
     * @param startPos
     * @param endPos
     * @return true if the node is already part of the graph
     */
    bool nodeExists(size_t startPos, size_t endPos) const;

    /**
     * get pointer to node in splice graph with same start and end position
     * @param startPos
     * @param endPos
     * @return pointer to node, or null if no such node in graph
     */
    SGNode* getNode(size_t startPos, size_t endPos) const;

    /**
     * get pointer to node in splice graph with nodeName
     * @param nodeName startPos-stopPos
     * @return pointer to node, or null if no such node in graph
     */
    SGNode* getNode(const string& nodeName) const;

    /**
     * get pointer to node in splice graph to which nodeName (partially) belongs
     * @param nodeName
     * @return pointer to node, or null if no such node exists
     */
    SGNode* getNodePartial(const string& nodeName) const;

    /**
     * get pointer to edge in splice graph with same small and large position
     * @param smallPos
     * @param largePos
     * @return pointer to edge, or null if no such edge in graph
     */
    SGEdge* getEdge(size_t smallPos, size_t largePos) const;

    /**
     * get pointer to edge in splice graph with edgename
     * @param edgeName smallPos_largePos
     * @return
     */
    SGEdge* getEdge(const string& edgeName) const;

    /**
     * Remove graph and corresponding nodes and edges
     */
    void deleteGraph();

    /**
     * Split node into 2 distinct nodes and adapt corresponding edges
     * The incoming edges now end in the new first node
     * The outgoing edges now start from the new second node
     * @param node node that needs to be split
     * @param splitPos startposition of second partial node
     * @return vector with the 2 new partial nodes
     */
    vector<SGNode*> splitNode(SGNode* node, size_t splitPos);

    /**
     * Add new junction to graph
     * @param junctionName using template chrName_smallPos_largePos_strand
     * @return true if new junction created, false if it already existed
     */
    bool addJunction(const string& junctionName);

    /**
     * get the node that is attached to a junction starting at position pos (+ create node if it does not yet exist)
     * @param pos smallest position of junction
     * @return
     */
    SGNode* getSmallNode(size_t pos);

    /**
    * get the node that is attached to a junction ending at position pos (+ create node if it does not yet exist)
    * @param pos largest position of junction
    * @return
    */
    SGNode* getLargeNode(size_t pos);

    /**
     * remove edge from graph
     * @param edgeToRemove
     */
    void removeEdge(SGEdge* edgeToRemove);

    /**
     * remove node from graph
     * @param node
     */
    void removeNode(SGNode* nodeToRemove);

    /**
     * merge consecutive nodes into 1 node if possible
     * @return the number of times a merge operation is performed
     */
    size_t mergeNodes();

    /**
     * Check if the graph can be alternatively spliced
     * @return
     */
    bool isSpliced();

    /**
     * get all nodes in interval [start, end]
     * @param start
     * @param end
     * @return sorted list of nodes in interval (only nodes that are completely inside interval)
     */
    vector<SGNode*> getNodesInInterval(size_t start, size_t end);

    /**
     * Get node in this graph that overlaps with the given consecutive edge
     * @param edge
     * @return
     */
    SGNode* getOverlappingNode(SGEdge* edge) const;

    /**
     * get pointer to edge between startNode and endNode
     * @param startNode
     * @param endNode
     * @return edge between start and endNode, NULL if edge does not exist
     */
    SGEdge* getExistingEdge(SGNode* startNode, SGNode* endNode) const;

    /**
     * merge graphs such that the resulting graph contains all the exons and junctions from both graphs
     *  Only apply to graphs representing the same gene
     * @param graph1
     * @param graph2
     * @result merged graph containing all the exons and junctions from both graphs
     */
    friend SpliceGraph mergeGraphs(const SpliceGraph& graph1, const SpliceGraph& graph2);

    /**
     * get mapping of exonnames and junctionnames of larger graph (this) to subgraph (graph2)
     * @param graph2
     * @return map of exon and junctionnames: if graph2 does not contain an exon/junction, it gets value ""
     */
    map<string, string> getVarMap(SpliceGraph& graph2) const;

    /**
     * write graph to file
     * @param fileName
     */
    void writeToFile(const string& fileName);

    /**
     * Fill observedNodes with the nodes in the graph that are observed in the depth file
     * @param exonDepthFile file containing exon depths for this gene
     * @param observedNodes reference to list of observed nodes
     */
     void getObservedNodes(const string& exonDepthFile, set<SGNode*>& observedNodes) const;

    /**
     * Fill observedEdges with the edges in the graph that are observed in the count file
     * @param junctionCountFile file containing junctin counts for this gene
     * @param observedEdges reference to list of observed edges
     * @param observedNodes reference to list of (previously computed) observed nodes
     * @param minJunctionCount minJunctionCount minimum depth for an observed junction with unobserved start/end exon(s) to be considered expressed
     * @param minJunctionCells minimum number of cells where an observed junction with unobserved start/end exon(s) is expressed to be included in the splicegraph
     */
    void getObservedEdges(const string& junctionCountFile, set<SGEdge*>& observedEdges, set<SGNode*>& observedNodes,
                          size_t  minJunctionCount, size_t minJunctionCells) const;


    /**
     * Make sure observedNodes and observedEdges make up a completely connected graph
     * @param observedNodes reference to list of observed nodes
     * @param observedEdges reference to list of observed edges
     */
    void makeConnected(set<SGNode*>& observedNodes, set<SGEdge*>& observedEdges);

    /**
     * Get tab separated varNames
     * @param orderedNodes contains list of ordered nodes as result
     * @param orderedEdges contains list of ordered edges as result
     * @return tab separated (ordered) nodeNames and edgeNames
     */
    string getVarNames(vector<SGNode*>& orderedNodes, vector<SGEdge*>& orderedEdges);

    /**
     * get number of paths from src to node
     * @param node
     * @return
     */
    uint1024_t getNumberOfPathsFromSrc(SGNode* node);


};

class VarGraph{
private:
    string geneID;
    size_t longestTranscriptLen;
    vector<Node*> nodes;
    vector<Edge*> edges;
    Node* srcNode;
    Node* sinkNode;

    size_t numVars;
    size_t numNonArtificialEdges;
    vector<string> varNames;

    map<SGNode*, Node*> nodeMap;
    map<SGEdge*, Edge*> edgeMap;

    map<Node*, vector<Var*>> spanningEdges;

public:
    explicit VarGraph(SpliceGraph* spliceGraph);

    void deleteGraph();

    size_t getNumVars() const { return numVars; };

    string getGeneID() const { return geneID; };

    size_t getLongestTranscriptLen() const { return longestTranscriptLen; };

    size_t getNumNonArtificialEdges() const { return numNonArtificialEdges; };

    Node* getSrcNode() const { return srcNode; };

    Node* getSinkNode() const { return sinkNode; };

    Var* getVar(SGNode* node) { return nodeMap[node]; };

    Var* getVar(SGEdge* edge);

    Var* getVar(const string& varName);

    size_t getVarID(SGNode* node) { return nodeMap[node]->getVarID(); };

    vector<Node*> getNodes() const{ return nodes; };

    vector<Edge*> getEdges() const{ return edges; };

    vector<Var*> getVars(){
        vector<Var*> vars(nodes.begin(), nodes.end());
        vector<Var*> edgeVars(edges.begin(), edges.end());
        vars.insert(vars.end(), edgeVars.begin(), edgeVars.end());
        return vars;
    };

    size_t getNumNodes() { return nodes.size(); };

    vector<size_t> getNonArtificialVarIDs();

    string getVarName(Var* var) { return varNames[var->getVarID()]; };

    string getVarName(size_t varID) { return varNames[varID]; };

    /**
     * get edges that span the node (start before node & end after node)
     * @param node
     * @return
     */
    vector<Var*> getSpanningEdges(Node* node);

    /**
     * get total number of paths from source to sink
     * @return
     */
    uint1024_t getTotalNumberOfPaths();

    /**
     * get number of paths from src to given node
     * @param node
     * @param numPathsFromSrc vector with already computed number of paths (contains 0 if not yet computed)
     * @return
     */
    uint1024_t getNumberOfPathsFromSrc(Node* node, vector<uint1024_t> &numPathsFromSrc);

};



#endif //ELLIPSIS_SPLICEGRAPH_H
