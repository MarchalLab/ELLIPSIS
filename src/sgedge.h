
#ifndef ELLIPSIS_SGEDGE_H
#define ELLIPSIS_SGEDGE_H

#include <vector>
#include <map>

using namespace std;

// ============================================================================
// SPLICE GRAPH EDGE CLASS
// ============================================================================


class SGNode; //forward declaration for pointers

class SGEdge {

private:
    SGNode* startNode;
    SGNode* endNode;
    string strand;

public:
    SGEdge(SGNode* startNode, SGNode* endNode, const string& strand);

    SGNode* getStartNode();
    SGNode* getEndNode();
    string getStrand();
    void setStartNode(SGNode* newStartNode);
    void setEndNode(SGNode* newEndNode);

    /**
     * get the smallest position of the corresponding junction
     * @return
     */
    size_t getEdgeSmallPos() const;

    /**
     * get the largest position of the corresponding junction
     * @return
     */
     size_t getEdgeLargePos() const;

    /**
     * get length of edge (= intron length)
     * @return
     */
    size_t getLength();

    /**
     * check if edge is an artificial edge (coming from source, going to sink or connecting 2 consecutive nodes)
     * @return true if edge is artificial
     */
    bool isArtificial();

    /**
     * check if ecge connects 2 consecutive nodes
     * @return true if edge connects consecutive nodes
     */
    bool isConsecutive();

    /**
     * check if edge is connected to source or sink node
     * @return true if edge is connected to source/sink node
     */
    bool isSourceSink() const;

    /**
     * get edge name
     * @return format startNode.nodeName -> endNode.nodeName
     */
    string getName() const;

    /**
     * get edge name
     * @return format smallPos_largePos
     */
    string getPosName() const;

    /**
     * Compare 2 edges
     * @param edge2
     * @return True if same start and endNode
     */
    bool operator==(const SGEdge& edge2);

    /**
     * Operator < overloading
     * @return true or false
     */
    bool operator<(const SGEdge& edge2) const {
        if (getEdgeSmallPos() == edge2.getEdgeSmallPos())
            return getEdgeLargePos() < edge2.getEdgeLargePos();
        return getEdgeSmallPos() < edge2.getEdgeSmallPos();
    }

    /**
     * operator << overloading
     * @param out output stream
     * @param e edge to display
     * @return output stream
     */
    friend std::ostream &operator<<(std::ostream &out, const SGEdge& e);

};

/**
 * Compare edges based on startNodes
 * @param edge1 Pointer to 1st edge
 * @param edge2 Pointer to 2nd edge
 * @return True if edge1->startNode < edge2->startNode
 */
bool compareEdges(SGEdge* edge1, SGEdge* edge2);

#endif //ELLIPSIS_SGEDGE_H
