
#ifndef ELLIPSIS_NNLSSOLVER_H
#define ELLIPSIS_NNLSSOLVER_H

#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Vector<double, Eigen::Dynamic> VectorXd;

/**
 * Class solving Non-Negative Least squares Ax=b
 */
class NNLSSolver {

private:
    MatrixXd AtA;
    MatrixXd At;
    size_t maxIterations;   //similar to Eigen/Unsupported/NNLS

    double maxVal;
    size_t iterations;
    size_t numEq;

public:
    NNLSSolver() = default;

    explicit NNLSSolver(MatrixXd A);

    NNLSSolver(MatrixXd A, double maxVal);

    /**
     * Solve Ax = b using non-negative least squares
     * @param b
     * @return x
     */
    VectorXd solve(VectorXd &b);

    /**
     * Solve Ax = b using non-negative least squares, starting at warmStart
     * @param b
     * @param warmStart
     * @return x
     */
    VectorXd solve(VectorXd &b, VectorXd &warmStart);

    size_t getIterations() const {
        return iterations;
    }

    size_t getNumEq() const{
        return numEq;
    }

};



#endif //ELLIPSIS_NNLSSOLVER_H
