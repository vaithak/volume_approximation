//
// Created by panagiotis on 25/6/2019.
//

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "Eigen"
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <vector>
#include <Eigen/Eigen>
#include <limits>

const double ZERO = 0.0000000001;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
typedef Eigen::SparseMatrix<double> SpMat;

/**
 * A linear matrix inequality A_0 + sum(x * F_i)
 */
class LMI {
    MT A0;

    // the matrices A_i, i>0
    std::vector<MT> matrices;

    typedef std::vector<MT>::iterator Iter;

public:
    LMI() {};

    LMI(std::vector<MT>& matrices) {
        this->A0 = matrices[0];

        for (int i=1 ; i<matrices.size() ; i++)
            this->matrices.push_back(matrices[i]);
    }

    LMI(const LMI& lmi) {
        this->A0 = lmi.A0;
        this->matrices = lmi.matrices;
    }

    /**
     * Evaluate the lmi for vector x
     *
     * @param x
     * @return
     */
    MT evaluate(const VT& x) {
       MT res = A0;
       int i = 0;

       for (Iter iter=matrices.begin() ; iter!=matrices.end() ; iter++, i++)
           res += x(i) * (*iter);

       return res;
    }

    bool isSingular(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() == 0)
                return true;

        return false;
    }

    bool isSingular(VT& x, double approx) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (abs(eivals(i).real()) <= abs(approx))
                return true;

        return false;
    }

    bool isNegativeDefinite(const VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() >= 0)
                return false;

        return true;
    }

    bool isNegativeSemidefinite(const VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++) {
            std::cout << eivals(i).real() << "\n";
            if (eivals(i).real() > 0)
                return false;
        }

        return true;
    }

    bool isPositiveSemidefinite(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() < 0)
                return false;

        return true;
    }

    bool isPositiveDefinite(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() <= 0)
                return false;

        return true;
    }


    bool isPositiveSemidefinite(VT& x, MT& mt, VT& minEigenVector, double& minEigenvalue) {
        bool res = true;

        mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();
        minEigenVector.resize(solver.eigenvectors().col(0).rows());

        for (int i=0 ; i<minEigenVector.rows() ; i++) {
            minEigenVector(i) = solver.eigenvectors().col(0)(i).real();

        }

        double min = eivals(0).real() / minEigenVector.norm();
        int minIndex = 0;

        for (int i = 0; i < eivals.rows(); i++) {
            for (int j=0 ; j<minEigenVector.rows() ; j++) {
                minEigenVector(j) = solver.eigenvectors().col(i)(j).real();
            }

            if (eivals(i).real() < 0)
                res = false;

            if (eivals(i).real() / minEigenVector.norm() < min) {
                min = eivals(i).real() / minEigenVector.norm();
                minIndex = i;
            }
        }


        minEigenvalue = eivals(minIndex).real()  / minEigenVector.norm() ;

        for (int i=0 ; i<minEigenVector.rows() ; i++) {
            minEigenVector(i) = solver.eigenvectors().col(minIndex)(i).real();

        }
        minEigenVector.normalize();

        return res;
    }


    int getDim() const {
        return matrices.size();
    }

    int getMatricesDim() const {
        return A0.rows();
    }

    const std::vector<MT>& getMatrices() const {
        return matrices;
    }

    /**
     * Evaluate the lmi for vector x without taking int account matrix A0
     *
     * @param x
     * @return
     */
    MT evaluateWithoutA0(const VT& x) {
        long dim = A0.rows();
        MT res;
        res.setZero(dim, dim);
        int i = 0;

        for (Iter iter=matrices.begin() ; iter!=matrices.end() ; iter++, i++)
            res += x(i) * (*iter);

        return res;
    }

    const MT& getAi(int i) const {
        return matrices[i];
    }

    const MT& getA0() const {
        return A0;
    }

    void setA0(MT& A0) {
        this->A0 = A0;
    }

    void addMatrix(MT& matrix) {
        matrices.push_back(matrix);
    }

    void print() {
        std::cout << "F0" << "\n" << A0 << "\n";
        int i = 1;

        for (Iter iter=matrices.begin() ; iter!=matrices.end() ; iter++, i++) {
            std::cout << "F" << i << "\n";
            std::cout << *iter << "\n";
        }
    }

    void changeInequalityDirection() {
        A0 = -1 * A0;

        for (int i=0 ; i<matrices.size() ; i++)
            matrices[i] = -1 * matrices[i];
    }

};


class Spectrahedron {
    /**
     * The collection of matrices that constitute the linear matrix
     * inequality describing the spectrahedron
     */
    LMI lmi;

    double maxDouble = std::numeric_limits<double>::max();
    double minDouble = std::numeric_limits<double>::lowest();

public:

    Spectrahedron() {}

    Spectrahedron(const Spectrahedron& spectrahedron) {
        LMI lmi;
        this->lmi = LMI(spectrahedron.getLMI());
    }

    Spectrahedron(LMI& lmi) {
        this->lmi = lmi;
    }

    const LMI& getLMI() const {
        return lmi;
    }



    /**
     * Compute the intersection of a 1D line and the spectrahedron by finding the
     * generalized eigenvalues of (LMI(x), A0 - LMI(direction))
     *
     * @param position
     * @param direction
     * @return (minimum positive eigenvalue, maximum negative eigenvalue)
     */
    std::pair<double, double> boundaryOracle(const VT& position, const VT& direction) {
        MT A = lmi.evaluate(position);
        MT B = -lmi.evaluateWithoutA0(direction);


        Eigen::GeneralizedEigenSolver<MT> ges(A,B);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();


        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;


        for (int i=0 ; i<alphas.rows() ; i++) {
            if (betas(i) == 0) //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative =lambda;
        }

        // for numerical stability
//        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
        if (lambdaMinPositive ==  maxDouble) lambdaMinPositive = 0; //TODO b must be too small..
        if (lambdaMaxNegative == minDouble) lambdaMaxNegative = 0;

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    /**
    * Compute the intersection of a 1D line and the spectrahedron by finding the
    * generalized eigenvalues of (LMI(x), A0 - LMI(direction))
    *
    * Take also in account the halfspace ax<=b
     *
    * @param position
    * @param direction
    * @return (minimum positive eigenvalue, maximum negative eigenvalue)
    */
    std::pair<double, double> boundaryOracleEfficient(const VT& position, const VT& direction, const VT& a, double b) {
        MT A = lmi.evaluate(position);
        MT B = -lmi.evaluateWithoutA0(direction);

        // Construct matrix operation object using the wrapper classes
        Spectra::DenseSymMatProd<double> op(B);
        Spectra::DenseCholesky<double>  Bop(-A);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, 3);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();
        // Retrieve results
        Eigen::VectorXd evalues;
//        Eigen::MatrixXd evecs;
        if(geigs.info() == Spectra::SUCCESSFUL)
        {
            evalues = geigs.eigenvalues();
//            evecs = geigs.eigenvectors();
        }

        double lambdaMaxNegative;
        double lambdaMinPositive;

        if (nconv == 2) {
            lambdaMaxNegative = - 1/evalues(0);
            lambdaMinPositive = -1/evalues(1);
        } else {
            lambdaMaxNegative = lambdaMinPositive = 0;
        }

        // check the cutting plane
        double lambda = (b - a.dot(position)) / a.dot(direction);
        if (lambda > 0 && lambda < lambdaMinPositive)
            lambdaMinPositive = lambda;
        if (lambda < 0 && lambda > lambdaMaxNegative)
            lambdaMaxNegative = lambda;

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    /**
    * Compute the intersection of a 1D line and the spectrahedron by finding the
    * generalized eigenvalues of (LMI(x), A0 - LMI(direction))
    *
    * Take also in account the halfspace ax<=b
     *
    * @param position
    * @param direction
    * @return (minimum positive eigenvalue, maximum negative eigenvalue)
    */
    std::pair<double, double> boundaryOracle(const VT& position, const VT& direction, const VT& a, double b) {
        MT A = lmi.evaluate(position);
        MT B = -lmi.evaluateWithoutA0(direction);

        Eigen::GeneralizedEigenSolver<MT> ges(A, B);

        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();

        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;

        for (int i=0 ; i<alphas.rows() ; i++) {

            if (betas(i) == 0)  //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative =lambda;
        }


        // for numerical stability
//        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
        if (lambdaMinPositive ==  maxDouble) lambdaMinPositive = 0; //TODO b must be too small..
        if (lambdaMaxNegative == minDouble) lambdaMaxNegative = 0;


        // check the cutting plane
        double lambda = (b - a.dot(position)) / a.dot(direction);
        if (lambda > 0 && lambda < lambdaMinPositive)
            lambdaMinPositive = lambda;
        if (lambda < 0 && lambda > lambdaMaxNegative)
            lambdaMaxNegative = lambda;

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    void print() {
        this->lmi.print();
    }

    void changeInequalityDirection() {
        lmi.changeInequalityDirection();
    }

    bool isSingular(VT& x) {
        return lmi.isSingular(x);
    }


    bool isSingular(VT& x, double approx) {
        return lmi.isSingular(x, approx);
    }

    std::pair<double, bool> boundaryOraclePositive(const VT& position, const VT& direction, const VT& a, double b, MT& LMIatP, MT& B, VT& genEivector, bool first = true) {
        if (first)
            LMIatP = lmi.evaluate(position);

//        if (!lmi.isNegativeSemidefinite(position)) throw "out\n";

        B = -lmi.evaluateWithoutA0(direction);

        // Construct matrix operation object using the wrapper classes
        Spectra::DenseSymMatProd<double> op(B);
        Spectra::DenseCholesky<double>  Bop(-LMIatP);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, 3);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();
        // Retrieve results
        Eigen::VectorXd evalues;
        Eigen::MatrixXd evecs;
        if(geigs.info() == Spectra::SUCCESSFUL)
        {
            evalues = geigs.eigenvalues();
            evecs = geigs.eigenvectors();
        }

        double lambdaMinPositive;

        if (nconv == 2) {
            lambdaMinPositive = -1/evalues(1);
            genEivector = evecs.col(1);
        } else {
            lambdaMinPositive = 0;
        }

        // check the cutting plane
        bool hitCuttingPlane = false;
        double lambda = (b - a.dot(position)) / a.dot(direction);

        if (lambda > 0 && lambda < lambdaMinPositive) {
            lambdaMinPositive = lambda;
            hitCuttingPlane = true;
        }
        B = -B;

        return {lambdaMinPositive, hitCuttingPlane};
    }

    double boundaryOracle_Boltzmann_HMC(const VT& position, const VT& direction, const VT& objectiveFunction, const double& temp,
            MT& B0, MT& B1, MT& B2, VT& genEivector, bool first = true) {

        unsigned int matrixDim= lmi.getMatricesDim();
        if (!lmi.isNegativeSemidefinite(position)) throw "out\n";
        lmi.print();
            std::cout << objectiveFunction << "\n";
//        if (first) {
            B0 = lmi.evaluate(position);
            B2 = lmi.evaluateWithoutA0(objectiveFunction);
//        }

        B1 = lmi.evaluateWithoutA0(direction);
        MT B2temp = B2 / (-2*temp);

        // create pencil matrix
        MT AA(2*matrixDim, 2*matrixDim);
        MT BB(2*matrixDim, 2*matrixDim);

        BB.block(matrixDim, matrixDim, matrixDim, matrixDim) = -1*B0;
        BB.block(0, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
        BB.block(matrixDim, 0, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
        BB.block(0, 0, matrixDim, matrixDim) = B2temp;

        AA.block(0, matrixDim, matrixDim, matrixDim) = B0;
        AA.block(0, 0, matrixDim, matrixDim) = B1;
        AA.block(matrixDim, 0, matrixDim, matrixDim) = B0;
        AA.block(matrixDim, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);


//        SpMat A =AA.sparseView();
//        SpMat B=BB.sparseView();
    std::cout << AA << "\n\n\n";
        std::cout << BB << "\n";

        // Construct matrix operation object using the wrapper classes
//        Spectra::SparseSymMatProd<double> op(A);
        Spectra::DenseSymMatProd<double> op(AA);
//        Spectra::SparseRegularInverse<double> Bop(B);

        Spectra::DenseCholesky<double>  Bop(-BB);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
//        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
//                geigs(&op, &Bop, 2, 4);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
//        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::SparseSymMatProd<double>, Spectra::SparseRegularInverse<double>, Spectra::GEIGS_REGULAR_INVERSE>
//                geigs(&op, &Bop, 1, 3);

        // Initialize and compute
//        geigs.init();
//        int nconv = geigs.compute();
//         Retrieve results
//        Eigen::VectorXd evalues;
//        Eigen::MatrixXd evecs;
//        if(geigs.info() == Spectra::SUCCESSFUL)
//        {
//            evalues = geigs.eigenvalues();
//            evecs = geigs.eigenvectors();
//        }
//
        double lambdaMinPositive = 0;

//        if (nconv == 2) {
//            lambdaMinPositive = evalues(0);
//            genEivector = evecs.col(0).segment(matrixDim, matrixDim);
//        } else {
//            lambdaMinPositive = 0;
//        }

//        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(AA, -1*BB);
//
//        bool _first = true;
//        int index;
//        const VT& eivals = es.eigenvalues();
//
//        for (int i=0 ; i<eivals.rows() ; i++) {
//            if (eivals(i) > 0 && (_first || eivals(i) < lambdaMinPositive)) {
//                lambdaMinPositive = eivals(i);
//                index = i;
//                _first = false;
//            }
//        }
//
//        MT vecs = es.eigenvectors();
//        genEivector = vecs.col(index).segment(matrixDim, matrixDim);

        Eigen::GeneralizedEigenSolver<MT> ges(AA, -BB);

        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();

        double lambdaMaxNegative = minDouble;
        lambdaMinPositive = maxDouble;
        int index = 0;

        for (int i=0 ; i<alphas.rows() ; i++) {

            if (betas(i) == 0 || alphas(i).imag() != 0)  //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive) {
                lambdaMinPositive = lambda;
                index = i;
            }
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative =lambda;
        }
        std::cout << lambdaMinPositive << "eval\n";

        Eigen::GeneralizedEigenSolver<MT>::EigenvectorsType eivecs = ges.eigenvectors();
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivec = eivecs.col(index);

        genEivector.setZero(matrixDim);

        for (int i = 0 ; i<matrixDim ; i++)
            genEivector(i) = eivec(matrixDim + i).real();

        std::cout << eivecs.col(index)<<"evec\n";
        return lambdaMinPositive;
    }

    void compute_reflection(const VT& genEivector, VT& direction, MT& C) {
        VT gradient;
        std::vector<MT> matrices = lmi.getMatrices();
        int dim = matrices.size();
        gradient = VT::Zero(dim);

        for (int i=0 ; i<dim ; i++) {
            gradient(i) = genEivector.dot((-matrices[i])* genEivector);
        }

//        Eigen::SelfAdjointEigenSolver<MT> es;
//        es.compute(C);
//        double product = 1;
//        VT evecs = es.eigenvalues();
//        for (int i=0 ; i<evecs.rows() ; i++) {
//            if (evecs(i)!= 0)
//                product *= evecs(i);
//        }

//        Eigen::FullPivLU<MT> lu_decomp(C);
//        auto rank = lu_decomp.rank();

//        std::cout << rank << "\n";
//        gradient *= product;

        gradient = -gradient;
        gradient.normalize();

        gradient = ((-2.0 * direction.dot(gradient)) * gradient);
        direction += gradient;
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
