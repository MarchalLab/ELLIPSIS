
#ifndef ELLIPSIS_VAR_H
#define ELLIPSIS_VAR_H

#include <vector>
#include <stdexcept>

using namespace std;

class Var {
private:
    size_t varID;
    bool isSrcSink;

public:
    vector<Var*> startVars;
    vector<Var*> endVars;

    Var(size_t varID, bool isSrcSink) : varID(varID), isSrcSink(isSrcSink) { };

    size_t getVarID() const { return varID; };

    bool isSrcSinkVar() const {return isSrcSink; };

    vector<Var*> getStartVars() const { return startVars; };

    vector<Var*> getEndVars() const { return endVars; };

    virtual bool isNode() const { return false; }; //entirely virtual function (should never be executed)

    virtual size_t getLength() const { return 0; }; //entirely virtual function (should never be executed)

};

class Node : public Var{
private:
    size_t len;

public:

    Node(size_t varID, size_t len, bool isSrcSink): Var(varID, isSrcSink), len(len) { };

    /**
     * Add incoming edge to this node
     * @param edge
     */
    void addInEdge(Var* edge);

    /**
     * Add outgoing edge to this node
     * @param edge
     */
    void addOutEdge(Var* edge);

    bool isNode() const override { return true; };

    size_t getLength() const override { return len; };

};


class Edge : public Var{
private:
    bool isConsecutive;

public:

    Edge(size_t varID, Node* startNode, Node* endNode, bool isSrcSink, bool isConsecutive);

    bool isNode() const override { return false; };

    size_t getLength() const override{ return 0; };

    /**
     * Check if edge connects 2 consecutive nodes
     * @return
     */
    bool isConsecutiveEdge() const { return isConsecutive; };

};

#endif //ELLIPSIS_VAR_H
