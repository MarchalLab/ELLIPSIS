
#include "PSIsolver.h"
#include "util.h"
#include <iostream>
#include <cassert>



PSIsolver::PSIsolver(VarGraph *graph, ObservedDepth *observedDepth, map<size_t, vector<size_t>> &neighbors,
                     int w_obs, int w_srcSink, int w_consv, int w_clust) :
    graph(graph), observedDepth(observedDepth), neighbors(neighbors), w_obs(w_obs), w_srcSink(w_srcSink),
    w_consv(w_consv), w_clust(w_clust), firstIter(true) {

    size_t nVar = graph->getNumVars();

    //already initialize fixed equations
    //get equations making sure src/sink node have PSI = 100%
    getSrcSinkEq();

    //get conservation of flow equations
    getConsvEq();

    //initialize X_obs (fixed values) & initialize y_obs (size)
    vector<size_t> nonArtificialVarIDs = graph->getNonArtificialVarIDs();
    X_obs = MatrixXd::Zero(int(nonArtificialVarIDs.size()), int(graph->getNumVars()));
    y_obs = VectorXd::Zero(int(nonArtificialVarIDs.size()));
    for (int i = 0; i < nonArtificialVarIDs.size(); i++){
        X_obs(i, int(nonArtificialVarIDs[i])) = w_obs;
    }

    //get solver if no neighborhood is observed
    //concatenate X
    size_t numEq = y_srcSink.size() + y_consv.size() + y_obs.size();
    MatrixXd X_noNB = MatrixXd(numEq, nVar);
    X_noNB << X_srcSink,
            X_consv,
            X_obs;
    //initialize Non-Negative Least Squares solver
    solver_noNB = NNLSSolver(X_noNB, 1);

    //initialize X_clust (fixed values) & initialize y_clust (size)
    X_clust = MatrixXd::Identity(int(nVar), int(nVar)) * w_clust;
    y_clust = VectorXd(nVar);

    //get solver if neighborhood is observed
    numEq += y_clust.size();
    //concatenate X
    MatrixXd X_NB = MatrixXd(numEq, nVar);
    X_NB << X_srcSink,
            X_consv,
            X_obs,
            X_clust;
    //NNLS solver
    solver_NB = NNLSSolver(X_NB, 1);

    oldPSI = vector<VectorXd>(observedDepth->getNumCells());

}

void PSIsolver::getSrcSinkEq(){

    //number of variables
    size_t nVar = graph->getNumVars();

    //number of equations
    int nEq = 2;

    //set equations for PSI(src) = 1 and PSI(sink) = 1
    X_srcSink = MatrixXd::Zero(nEq, int(nVar));
    y_srcSink = VectorXd ::Constant(nEq, w_srcSink);

    X_srcSink(0, int(graph->getSrcNode()->getVarID())) = w_srcSink;
    X_srcSink(1, int(graph->getSinkNode()->getVarID())) = w_srcSink;
}

void PSIsolver::getConsvEq(){

    //number of variables
    size_t nVar = graph->getNumVars();

    //2 equations per node (incoming + outgoing), except for src/sink
    int nEq = (int(graph->getNumNodes()) - 1) * 2;

    //initialize equations
    X_consv = MatrixXd::Zero(nEq, int(nVar));
    y_consv = VectorXd::Zero(nEq);

    //compute conservation of flow equations
    int eqCtr = 0;
    for (auto const &node : graph->getNodes()){
        int varID = int(node->getVarID());

        //add equations for incoming junctions
        vector<Var*> inEdges = node->getStartVars();
        if (! inEdges.empty()){
            X_consv(eqCtr, varID) = w_consv;
            for (auto const &edge : inEdges)
                X_consv(eqCtr,int(edge->getVarID())) = -w_consv;
            eqCtr++;
        }

        //add equations for outgoing junctions
        vector<Var*> outEdges = node->getEndVars();
        if (! outEdges.empty()){
            X_consv(eqCtr, varID) = w_consv;
            for (auto const &edge : outEdges)
                X_consv(eqCtr,int(edge->getVarID())) = -w_consv;
            eqCtr++;
        }

    }

}

void PSIsolver::getObsEq(int cellIdx){

    int eqCtr = 0;

    //add equation for each observed node
    double alphaExon = observedDepth->getAlphaExon(cellIdx);
    for (auto const &node : graph->getNodes()){
        //skip src/sink node
        if (node->isSrcSinkVar())
            continue;

        int varID = int(node->getVarID());
        assert(eqCtr >= 0 and eqCtr < y_obs.size());
        y_obs(eqCtr) = observedDepth->getExonDepth(cellIdx, varID) / alphaExon * w_obs ;

        eqCtr++;
    }

    //add equation for each observed edge
    double alphaJunction = observedDepth->getAlphaJunction(cellIdx);
    for (auto const &edge : graph->getEdges()){
        //skip src/sink edges
        if (edge->isSrcSinkVar())
            continue;

        int varID = int(edge->getVarID());
        assert(eqCtr >= 0 and eqCtr < y_obs.size());
        if (edge->isConsecutiveEdge())
            y_obs(eqCtr) = observedDepth->getJunctionDepth(cellIdx, varID) / alphaExon * w_obs ;
        else
            y_obs(eqCtr) = observedDepth->getJunctionDepth(cellIdx, varID) / alphaJunction * w_obs;

        eqCtr++;
    }

}

void PSIsolver::getClusterEq(int cellIdx){

    for (auto const &node : graph->getNodes()){
        int varID = int(node->getVarID());
        double sumPSI = 0;
        double sumAlpha = 0;

        for (auto const &neighbor : neighbors[cellIdx]){

            //cells with higher depth (= high alpha) have more accurate PSI estimates => weight using alpha
            assert(varID >= 0 and varID < oldPSI[neighbor].size());
            sumPSI += oldPSI[neighbor][varID] * observedDepth->getAlphaExon(neighbor);
            sumAlpha += observedDepth->getAlphaExon(neighbor);

        }

        assert(varID >= 0 and varID < y_clust.size());
        y_clust(varID) = sumPSI / sumAlpha * w_clust;
    }

    for (auto const &edge : graph->getEdges()){
        int varID = int(edge->getVarID());
        double sumPSI = 0;
        double sumAlpha = 0;

        for (auto const &neighbor : neighbors[cellIdx]){

            //cells with higher depth (= high alpha) have more accurate PSI estimates => weight using alpha
            assert(varID >= 0 and varID < oldPSI[neighbor].size());
            sumPSI += oldPSI[neighbor][varID] * observedDepth->getAlphaJunction(neighbor);
            sumAlpha += observedDepth->getAlphaJunction(neighbor);

        }

        assert(varID >= 0 and varID < y_clust.size());
        y_clust(varID) = sumPSI / sumAlpha * w_clust;
    }

}

bool PSIsolver::solveEqs() {

    //number of cells that reached convergence
    int numCellsConv = 0;

    for (auto const &it : neighbors){

        int cellIdx = int(it.first);

        //get depth equations
        getObsEq(cellIdx);

        VectorXd newPSI;

        if (! firstIter and ! neighbors[cellIdx].empty()){

            //get cluster equations
            getClusterEq(cellIdx);

            VectorXd y(solver_NB.getNumEq());
            y << y_srcSink,
                    y_consv,
                    y_obs,
                    y_clust;

            newPSI = solver_NB.solve(y, oldPSI[cellIdx]);
        } else{

            VectorXd y(solver_noNB.getNumEq());
            y << y_srcSink,
                    y_consv,
                    y_obs;


            if (firstIter)
                newPSI = solver_noNB.solve(y);
            else
                newPSI = solver_noNB.solve(y, oldPSI[cellIdx]);
        }


        assert(newPSI.size() == graph->getNumVars());
        for (int var = 0; var < newPSI.size(); var++){
            if (isnan(newPSI[var]))
                cout << "newPSI = NaN" << endl;
        }

        //compute average cross section PSI sum
        double sumCrossSection = 0;
        for (auto node: graph->getNodes()) {
            double sumPSI = 0;

            if (!node->isSrcSinkVar())
                sumPSI += newPSI(int(node->getVarID()));

            for (auto edge: graph->getSpanningEdges(node)) {
                sumPSI += newPSI(int(edge->getVarID()));
            }

            sumCrossSection += sumPSI;
        }
        double avCrossSection = sumCrossSection / double(graph->getNumNodes());

        //rescale PSI values to make av cross section = 1
        newPSI = newPSI / avCrossSection;

        //max PSI value is 1
        newPSI = newPSI.cwiseMin(1);

        //check convergence : cell is converged if PSI for each variable differs less than 0.01% with prev. iteration
        if ( !firstIter and ((oldPSI[cellIdx] - newPSI).cwiseAbs().array() <= 0.0001).all() )
            numCellsConv++;

        //update oldPSI
        assert(cellIdx >= 0 and cellIdx < oldPSI.size());
        oldPSI[cellIdx] = newPSI;

    }

    //oldPSI is now updated
    firstIter = false;

    //PSI for this gene is converged if all cells are converged
    return numCellsConv == observedDepth->getNumCells();
}

vector<VectorXd> PSIsolver::getPSI(){
    return oldPSI;
}

void PSIsolver::writeResults(const string& outFile){

    //open file
    ofstream OFS(outFile);


    //write header
    OFS << "cell";
    for (int varID = 0; varID < graph->getNumVars(); varID++)
        OFS << "\t" << graph->getVarName(varID);
    OFS << endl;

    //write PSI for each cell
    for (int cellIdx = 0; cellIdx < oldPSI.size(); cellIdx++ ){

        OFS << observedDepth->getCell(cellIdx) ;
        for (int varID = 0; varID < graph->getNumVars(); varID++)
            OFS << "\t" << round(100 * oldPSI[cellIdx](varID)) / 100;
        OFS << endl;

    }

    //close file
    OFS.close();

}
