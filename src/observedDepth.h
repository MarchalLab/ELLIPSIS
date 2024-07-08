
#ifndef ELLIPSIS_OBSERVEDDEPTH_H
#define ELLIPSIS_OBSERVEDDEPTH_H

#include <string>
#include "spliceGraph.h"
#include <Eigen/Dense>

typedef Eigen::Vector<double, Eigen::Dynamic> VectorXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

class ObservedDepth {
private:
    VarGraph* graph;

    vector<vector<double>> cellVarDepth;    //observed depth for each cell & variable

    vector<double> cellAlphaInit;           //initial alpha estimate
    vector<double> cellAlphaExon;           //alpha for exons for each cell
    vector<double> cellAlphaJunction;       //alpha for junctions for each cell

    vector<string> cells;

    void getExonDepth(SpliceGraph* spliceGraph, const string& exonDepthFile);
    void getJunctionDepth(SpliceGraph* spliceGraph, const string& junctionDepthFile);

public:
    ObservedDepth(SpliceGraph* spliceGraph, VarGraph* graph, const string& exonDepthFile, const string& junctionDepthFile);

    /**
     * initialize alpha (cell specific depth) using the exons with maximum depth
     */
    void initAlpha();

    /**
     * get overall cell specific depth value alpha (for exons)
     * @param cellIdx
     * @return
     */
    double getAlphaExon(size_t cellIdx);

    /**
     * get overall cell specific depth value alpha (for junctions)
     * @param cell
     * @return
     */
    double getAlphaJunction(size_t cellIdx);

    /**
     * Update values for cell specific depth alpha for each cell
     *  solving equation Xa = y using Non-Negative Least squares
     *  with X = vector containing the previously computed PSI values, and y the (scaled) observed depths
     * @param cellPSI
     */
    void updateAlpha(vector<VectorXd> &cellPSI);

    /**
     * get (scaled) depth for cell and varID
     * @param cellIdx
     * @param varID
     * @return observed depth
     */
    double getExonDepth(size_t cellIdx, size_t varID);

    /**
     * get (scaled) depth for cell and varID
     * @param cellIdx
     * @param varID
     * @return observed depth
     */
    double getJunctionDepth(size_t cellIdx, size_t varID);

    /**
     * Get cell name for cell index
     * @param cellIdx
     * @return
     */
    string getCell(size_t cellIdx);

    /**
     * Get number of cells
     * @return
     */
    size_t getNumCells();

    /**
     * Remove cells for which alpha < minDepth
     * @param minDepth depth limit for which cells need to be filtered out
     * @param lowQualFrac low quality read fractions per cell
     * @param maxLowQual maximum fraction of low quality read mappings
     * @param neighbors
     * @param outFile
     */
    map<size_t, vector<size_t>> filterCells(double minDepth, map<string, double> lowQualFrac, double maxLowQual,
                                      map<string, vector<string>>& neighbors, const string& outFile);

    /**
     * Write cell specific alpha to output stream
     * @param outFile
     */
    void writeAlpha(const string& outFile);

    /**
     * get cellNames
     * @return
     */
    vector<string> getCells(){ return cells; }

    /**
     * write number of (expressed) neighbors to file
     * @param outFile
     * @param neighbors
     */
    void writeNumNB(const string& outFile, const map<size_t, vector<size_t>>& neighbors);

};


#endif //ELLIPSIS_OBSERVEDDEPTH_H
