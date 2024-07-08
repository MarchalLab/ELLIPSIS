
#ifndef ELLIPSIS_PSISOLVER_H
#define ELLIPSIS_PSISOLVER_H

#include <string>
#include <map>
#include "observedDepth.h"
#include "observedDepth.h"
#include "NNLSSolver.h"
#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Vector<double, Eigen::Dynamic> VectorXd;

using namespace std;

class PSIsolver {
private:
    VarGraph* graph;
    ObservedDepth* observedDepth;
    map<size_t, vector<size_t>> &neighbors;
    bool firstIter;

    //equations
    MatrixXd X_srcSink;
    VectorXd y_srcSink;
    int w_srcSink;

    MatrixXd X_consv;
    VectorXd y_consv;
    int w_consv;

    MatrixXd X_obs;
    VectorXd y_obs;
    int w_obs;

    MatrixXd X_clust;
    VectorXd y_clust;
    int w_clust;

    //Non-Negative Least Squares solver
    NNLSSolver solver_NB;   //solver if there are neighbors observed
    NNLSSolver solver_noNB; //solver if there are neighbors observed

    vector<VectorXd> oldPSI;

    /**
     * get (weighted) equations source and sink node
     *      PSI(src) = 1
     *      PSI(sink) = 1
     */
    void getSrcSinkEq();

    /**
     * get (weighted) equations for conservation of flow
     *     PSI(node) - PSI(incoming edges) = 0
     *     PSI(node) - PSI(outgoing edges) = 0
     */
    void getConsvEq();

    /**
     * get (weighted) equations for observed depth
     *      PSI(var) * alpha = (scaled) obsDepth
     */
    void getObsEq(int cellIdx);

    /**
     * get (weighted) equations for intra-cell type similarity
     *      PSI(var) = avg(PSI(var) for cells in neighborhood)
     * @param cluster
     */
    void getClusterEq(int cellIdx);

public:
    PSIsolver(VarGraph* graph, ObservedDepth* observedDepth, map<size_t, vector<size_t>> &neighbors, int wObs,
              int wSrcSink, int wFlow, int wClust);

    /**
     * get solution for all equations using the Non-Negative Least Squares method
     *  Ax = y : get x such that square error (Ax - y)² is minimized, while making sure all x >= 0
     *  @return true if solution converged (less than 1% change in all cells)
     */
    bool solveEqs();

    /**
     * get computed PSI values
     * @return PSI values for each cell
     */
    vector<VectorXd> getPSI();

    /**
     * Write results to file
     * @param outFile
     */
    void writeResults(const string& outFile);

};


#endif //ELLIPSIS_PSISOLVER_H
