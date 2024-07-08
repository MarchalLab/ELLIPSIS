
#include "NNLSSolver.h"

NNLSSolver::NNLSSolver(MatrixXd A){
    At = A.transpose();
    AtA = At * A;

    maxIterations = 2 * A.cols();
    maxVal = std::numeric_limits<int>::max(); //no upper limit for solutions

    numEq = A.rows();
}

NNLSSolver::NNLSSolver(MatrixXd A, double maxVal) : maxVal(maxVal) {
    At = A.transpose();
    AtA = At * A;

    maxIterations = 2 * A.cols();

    numEq = A.rows();
}

VectorXd NNLSSolver::solve(VectorXd &b, VectorXd &warmStart){
    //code adapted from https://codereview.stackexchange.com/questions/277691/coordinate-descent-non-negative-least-squares-optimization

    VectorXd Atb = At * b;
    size_t maxit = maxIterations;

    VectorXd x = warmStart;
    for (int i = 0; i < Atb.size(); i++)
        Atb -= AtA.col(i) * x(i);

    double tol = 1e-8 * double(Atb.size());

    while(maxit-- > 0) {
        double tol_ = 0;
        for (int i = 0; i < Atb.size(); ++i) {
            double diff = Atb(i) / AtA(i, i);

            if (-diff > x(i)) { //constraint x >= 0
                if (x(i) != 0) {
                    Atb -= AtA.col(i) * -x(i);
                    tol_ = 1;
                    x(i) = 0;
                }
            }else if (maxVal - diff < x(i)){ //constraint x <= maxVal
                if (x(i) != maxVal){
                    Atb -= AtA.col(i) * (maxVal-x(i));
                    tol_ = std::abs((maxVal-x(i)) / maxVal);
                    x(i) = maxVal;
                }
            }else if (diff != 0) {
                x(i) += diff;
                Atb -= AtA.col(i) * diff;
                tol_ += std::abs(diff / x(i));
            }
        }
        if (tol_ < tol) break;
    }
    iterations = maxIterations - maxit;
    return x;
}

VectorXd NNLSSolver::solve(VectorXd &b) {
    //code adapted from https://codereview.stackexchange.com/questions/277691/coordinate-descent-non-negative-least-squares-optimization

    VectorXd Atb = At * b;
    size_t maxit = maxIterations;

    VectorXd x = VectorXd::Zero(Atb.size());

    double tol = 1e-8 * double(Atb.size());

    while(maxit-- > 0) {
        double tol_ = 0;
        for (int i = 0; i < Atb.size(); ++i) {
            double diff = Atb(i) / AtA(i, i);

            if (-diff > x(i)) { //constraint x >= 0
                if (x(i) != 0) {
                    Atb -= AtA.col(i) * -x(i);
                    tol_ = 1;
                    x(i) = 0;
                }
            }else if (maxVal - diff < x(i)){ //constraint x <= maxVal
                if (x(i) != maxVal){
                    Atb -= AtA.col(i) * (maxVal-x(i));
                    tol_ = std::abs((maxVal-x(i)) / maxVal);
                    x(i) = maxVal;
                }
            }else if (diff != 0) {
                x(i) += diff;
                Atb -= AtA.col(i) * diff;
                tol_ += std::abs(diff / x(i));
            }
        }
        if (tol_ < tol) break;
    }
    iterations = maxIterations - maxit;
    return x;
}
