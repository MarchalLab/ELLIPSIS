
#ifndef ELLIPSIS_SGNODE_H
#define ELLIPSIS_SGNODE_H

#include <vector>
#include <map>

using namespace std;

// ============================================================================
// SPLICE GRAPH NODE CLASS
// ============================================================================

class SGEdge; //forward declaration for in/outedges

class SGNode {

private:
    vector<SGEdge*> inEdges;
    vector<SGEdge*> outEdges;

protected:
    size_t startPos;
    size_t endPos;

public:
    SGNode(size_t startPos, size_t endPos);
    explicit SGNode(const string& nodeName);
    SGNode(const SGNode& oldNode); //copy constructor

    size_t getStartPos() const;
    size_t getEndPos() const;
    void setStartPos(size_t newStartPos);
    void setEndPos(size_t newEndPos);
    vector<SGEdge*> getInEdges();
    vector<SGEdge*> getOutEdges();


    /**
     * Compare 2 nodes according to startPos
     * @param node2
     * @return True if node1.startPos < node2.startPos
     */
    bool operator<(const SGNode& node2) const;

    /**
     * Compare 2 nodes
     * @param node2
     * @return True if same start and endPos
     */
    bool operator==(const SGNode& node2) const;

    /**
     * check if the positions of node2 overlap with this node
     * @param node2
     * @return
     */
    bool overlaps(SGNode* node2) const;

    /**
     * check if this node starts before node2
     * @param node2
     * @return
     */
    bool startsBefore(SGNode* node2) const;

    /**
     * Check if this node ends before node2
     * @param node2
     * @return
     */
    bool endsBefore(SGNode* node2) const;

    /**
     * Add new incoming edge
     * @param newEdge
     */
    void addInEdge(SGEdge* newEdge);

    /**
     * Add new outgoing edge
     * @param newEdge
     */
    void addOutEdge(SGEdge* newEdge);

    /**
     * remove edge from list of incoming edges
     * @param edgeToRemove
     */
    void removeInEdge(SGEdge* edgeToRemove);

    /**
     * remove edge from list of outgoing edges
     * @param edgeToRemove
     */
    void removeOutEdge(SGEdge* edgeToRemove);

    /**
     * Get the number of bases of this node
     * @return
     */
    virtual size_t getLength();

    /**
     * Check if node is artificial node (source or sink node)
     * @return true if artificial
     */
    virtual bool isArtificial();

    /**
     * check if node is source node
     * @return true if source node
     */
    virtual bool isSource();

    /**
     * check if node is source node
     * @return true if source node
     */
    virtual bool isSink();

    /**
     * get node name
     * @return format startPos-endPos
     */
    virtual string getName();

    /**
     * operator << overloading
     * @param out output stream
     * @param n node to display
     * @return output stream
     */
    friend std::ostream &operator<<(std::ostream &out, const SGNode &n);

};

class SourceSGNode: public SGNode{

public:

    /**
     * Constructor
     * @param strand + or - strand ('.' = undefined)
     */
    explicit SourceSGNode(const string& strand);

    /**
     * Get the number of bases of this node
     * @return
     */
    size_t getLength() override;

    /**
     * Check if node is artificial node (source or sink node)
     * @return true, because source node is artificial
     */
    bool isArtificial() override;

    /**
     * check if node is source node
     * @return true
     */
    bool isSource() override;

    /**
     * check if node is sink node
     * @return false
     */
    bool isSink() override;

    /**
     * get node name
     * @return src
     */
    string getName() override;
};

class SinkSGNode: public SGNode{

public:

    /**
     * Constructor
     * @param strand + or - strand ('.' = undefined)
     */
    explicit SinkSGNode(const string& strand);

    /**
     * Get the number of bases of this node
     * @return
     */
    size_t getLength() override;

    /**
     * Check if node is artificial node (source or sink node)
     * @return true, because sink node is artificial
     */
    bool isArtificial() override;

    /**
     * check if node is source node
     * @return false
     */
    bool isSource() override;

    /**
     * check if node is sink node
     * @return true
     */
    bool isSink() override;

    /**
     * get node name
     * @return sink
     */
    string getName() override;
};

/**
 * Compare nodes based on startpostion
 * @param node1 Pointer to 1st node
 * @param node2 Pointer to 2nd node
 * @return True if node1->startPos < node2->startPos
 */
bool compareNodes(SGNode* node1, SGNode* node2);

#endif //ELLIPSIS_SGNODE_H
