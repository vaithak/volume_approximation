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
#include <chrono>
#include <limits>

const double ZERO = 0.0000000001;
double ORACLE_TIME;
double REFLECTION_TIME;
int BOUNDARY_CALLS;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MT_ROWMAJOR;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
typedef Eigen::SparseMatrix<double> SpMat;

/**
 * A linear matrix inequality A_0 + sum(x * F_i)
 */
class LMI {
    MT A0;

    // the matrices A_i, i>0
    std::vector<MT> matrices;
    MT_ROWMAJOR vectorMatrix;
    typedef std::vector<MT>::iterator Iter;

public:
    LMI() {};

    LMI(std::vector<MT> &matrices) {
        this->A0 = matrices[0];
        for (int i = 1; i < matrices.size(); i++) {
            this->matrices.push_back(matrices[i]);
        }

        setVectorMatrix();
    }

    LMI(const LMI &lmi) {
        this->A0 = lmi.A0;
        this->matrices = lmi.matrices;
        setVectorMatrix();
    }

    /**
     * Evaluate the lmi for vector x
     *
     * @param x
     * @return
     */
    void evaluate(const VT &x, MT &res) {
        int dim = A0.rows();
        res = MT::Zero(dim, dim);
        res += A0;
        int i = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++)
            res.noalias() += x(i) * (*iter);
    }

    void evaluate_revised(const VT &x, MT &res) {
        res.setConstant(0);
        res += A0;
        int i = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++)
            res.noalias() += x(i) * (*iter);
    }


    bool isSingular(VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() == 0)
                return true;

        return false;
    }

    bool isSingular(VT &x, double approx) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (abs(eivals(i).real()) <= abs(approx))
                return true;

        return false;
    }

    bool isNegativeDefinite(const VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() >= 0)
                return false;

        return true;
    }

    bool isNegativeSemidefinite(const VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++) {
//            std::cout << eivals(i).real() << "\n";
            if (eivals(i).real() > 0)
                return false;
        }

        return true;
    }

    bool isPositiveSemidefinite(VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() < 0)
                return false;

        return true;
    }

    bool isPositiveDefinite(VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() <= 0)
                return false;

        return true;
    }


    bool isPositiveSemidefinite(VT &x, MT &mt, VT &minEigenVector, double &minEigenvalue) {
        bool res = true;

        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();
        minEigenVector.resize(solver.eigenvectors().col(0).rows());

        for (int i = 0; i < minEigenVector.rows(); i++) {
            minEigenVector(i) = solver.eigenvectors().col(0)(i).real();

        }

        double min = eivals(0).real() / minEigenVector.norm();
        int minIndex = 0;

        for (int i = 0; i < eivals.rows(); i++) {
            for (int j = 0; j < minEigenVector.rows(); j++) {
                minEigenVector(j) = solver.eigenvectors().col(i)(j).real();
            }

            if (eivals(i).real() < 0)
                res = false;

            if (eivals(i).real() / minEigenVector.norm() < min) {
                min = eivals(i).real() / minEigenVector.norm();
                minIndex = i;
            }
        }


        minEigenvalue = eivals(minIndex).real() / minEigenVector.norm();

        for (int i = 0; i < minEigenVector.rows(); i++) {
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

    const std::vector<MT> &getMatrices() const {
        return matrices;
    }

    /**
     * Evaluate the lmi for vector x without taking int account matrix A0
     *
     * @param x
     * @return
     */
    void evaluateWithoutA0(const VT &x, MT &res) {
        long dim = A0.rows();

        res = MT::Zero(dim, dim);

        int i = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++)
            res.noalias() += x(i) * (*iter);
    }

    /**
 * Evaluate the lmi for vector x without taking int account matrix A0
 *
 * @param x
 * @return
 */
    void evaluateWithoutA0_revised(const VT &x, MT &res) {

        long dim = A0.rows();

        VT a = vectorMatrix * x;

        double *data = res.data();
        double *v = a.data();

        int at = 0;

        // copy lower triangular
        for (int at_col = 0; at_col < dim; at_col++) {
            int col_offset = at_col * dim;
            double *target = data + col_offset + at_col;

            for (int at_row = at_col; at_row < dim; at_row++) {
                *(target++) = *(v++);
            }
        }

        v = a.data();

        // copy upper triangular
        for (int at_row = 0; at_row < dim; at_row++) {
            double *target = data + at_row + at_row * dim;

            for (int at_col = at_row; at_col < dim; at_col++) {
                *target = *(v++);
                target = target + dim;
            }
        }

    }

    void setVectorMatrix() {
        int dim = getMatricesDim();
        int newDim = dim * (dim + 1) / 2;

        vectorMatrix.resize(newDim, getDim());

        int atMatrix = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, atMatrix++) {
            int i = 0;

            for (int at_row = 0; at_row < dim; at_row++)
                for (int at_col = at_row; at_col < dim; at_col++) {
                    vectorMatrix(i++, atMatrix) = (*iter)(at_row, at_col);
                }

        }

    }

    void evaluateWithoutA0_revised2(const VT &x, MT &res) {
        int i = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++)
            res.noalias() += x(i) * (*iter);
    }

    void evaluateWithoutA0(const VT &x, Eigen::TriangularView<MT, Eigen::Upper> &res) {
        long dim = A0.rows();

        res = MT::Zero(dim, dim);
        int i = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++)
            res += x(i) * (*iter);
    }

    const MT &getAi(int i) const {
        return matrices[i];
    }

    const MT &getA0() const {
        return A0;
    }

    void setA0(MT &A0) {
        this->A0 = A0;
    }

    void addMatrix(MT &matrix) {
        matrices.push_back(matrix);
    }

    void print() {
        std::cout << "F0" << "\n" << A0 << "\n";
        int i = 1;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++) {
            std::cout << "F" << i << "\n";
            std::cout << *iter << "\n";
        }
    }

    void changeInequalityDirection() {
        A0 = -1 * A0;

        for (int i = 0; i < matrices.size(); i++)
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

    Spectrahedron(const Spectrahedron &spectrahedron) {
        LMI lmi;
        this->lmi = LMI(spectrahedron.getLMI());
    }

    Spectrahedron(LMI &lmi) {
        this->lmi = lmi;
    }

    const LMI &getLMI() const {
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
    std::pair<double, double> boundaryOracle(const VT &position, const VT &direction) {
        MT A;
        lmi.evaluate(position, A);
        MT B;
        lmi.evaluateWithoutA0(-1 * direction, B);


        Eigen::GeneralizedEigenSolver<MT> ges(A, B);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();

        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;


        for (int i = 0; i < alphas.rows(); i++) {
            if (betas(i) == 0) //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative = lambda;
        }

        // for numerical stability
//        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
        if (lambdaMinPositive == maxDouble) lambdaMinPositive = 0; //TODO b must be too small..
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
    std::pair<double, double>
    boundaryOracleEfficient(const VT &position, const VT &direction, const VT &a, double b) {
        MT A;
        lmi.evaluate(position, A);
        MT B;
        lmi.evaluateWithoutA0(-1 * direction, B);
        BOUNDARY_CALLS++;
        // Construct matrix operation object using the wrapper classes
        Spectra::DenseSymMatProd<double> op(B);
        Spectra::DenseCholesky<double> Bop(-A);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, 3);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();
        // Retrieve results
        Eigen::VectorXd evalues;
//        Eigen::MatrixXd evecs;
        if (geigs.info() == Spectra::SUCCESSFUL) {
            evalues = geigs.eigenvalues();
//            evecs = geigs.eigenvectors();
        }

        double lambdaMaxNegative;
        double lambdaMinPositive;

        if (nconv == 2) {
            lambdaMaxNegative = -1 / evalues(0);
            lambdaMinPositive = -1 / evalues(1);
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
    std::pair<double, double> boundaryOracle(const VT &position, const VT &direction, const VT &a, double b) {
        MT A;
        lmi.evaluate(position, A);
        MT B;
        lmi.evaluateWithoutA0(-1 * direction, B);

        Eigen::GeneralizedEigenSolver<MT> ges(A, B);
        BOUNDARY_CALLS++;
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();

        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;

        for (int i = 0; i < alphas.rows(); i++) {

            if (betas(i) == 0)  //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative = lambda;
        }


        // for numerical stability
//        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
        if (lambdaMinPositive == maxDouble) lambdaMinPositive = 0; //TODO b must be too small..
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

    bool isSingular(VT &x) {
        return lmi.isSingular(x);
    }


    bool isSingular(VT &x, double approx) {
        return lmi.isSingular(x, approx);
    }

    template<class Point>
    class BoundaryOracleBilliardSettings {
    public:
        MT B, LMIatP;
        Point genEigenvector;
        double max_segment_length;
        bool first; //true if first call of the boundary oracle

        BoundaryOracleBilliardSettings() {
            first = true;
            max_segment_length = 0;
        }

        void setMaxSegmentLength(double lambda) {
            if (lambda > 0 && lambda > max_segment_length)
                max_segment_length = lambda;
        }

        void resetMaxSegmentLength() {
            max_segment_length = 0;
        }

        double maxSegmentLength() {
            return max_segment_length;
        }
    };

    template<class Point>
    std::pair<double, bool>
    boundaryOracleBilliard(const VT &position, const VT &direction, const VT &a, const double &b,
                           BoundaryOracleBilliardSettings<Point> &settings) {

        if (settings.first) {
            int dim = lmi.getMatricesDim();
            settings.LMIatP.resize(dim, dim);
            lmi.evaluate_revised(position, settings.LMIatP);
            settings.B.resize(dim, dim);
        }

        BOUNDARY_CALLS++;
//        if (!lmi.isNegativeSemidefinite(position)) throw "out\n";

        lmi.evaluateWithoutA0_revised(direction, settings.B);

        auto t1 = std::chrono::steady_clock::now();

        // Construct matrix operation object using the wrapper classes
        Spectra::DenseSymMatProd<double> op(settings.B);
        Spectra::DenseCholesky<double> Bop(-settings.LMIatP);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 1, 15);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        auto t2 = std::chrono::steady_clock::now();

        ORACLE_TIME += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

        // Retrieve results
        Eigen::VectorXd evalues;
        Eigen::MatrixXd evecs;

        if (geigs.info() == Spectra::SUCCESSFUL) {
            evalues = geigs.eigenvalues();
            evecs = geigs.eigenvectors();
        }

        double lambdaMinPositive;

        if (nconv == 1) {
            lambdaMinPositive = 1 / evalues(0);
            settings.genEigenvector = Point(evecs.col(0));
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

        settings.setMaxSegmentLength(lambdaMinPositive);


        return {lambdaMinPositive, hitCuttingPlane};
    }


    template<class Point>
    class BoundaryOracleBoltzmannHMCSettings {
    public:
        MT B0, B1, B2;
        Point genEigenvector;
        Point gradient;
        bool first; //true if first call of the boundary oracle
        double epsilon; //when a point is this close to the boundary, consider it a boundary point
        bool isBoundaryPoint;

        BoundaryOracleBoltzmannHMCSettings() {
            first = true;
            epsilon = 0.0001;
        }
    };

    template<class Point>
    double
    boundaryOracle_Boltzmann_HMC(const Point &_position, const Point &_direction, const Point &_objectiveFunction,
                                 const double &temp, BoundaryOracleBoltzmannHMCSettings<Point> &settings) {

        const VT &position = _position.getCoefficients();
        const VT &direction = _direction.getCoefficients();
        const VT &objectiveFunction = _objectiveFunction.getCoefficients();

        unsigned int matrixDim = lmi.getMatricesDim();
        if (!lmi.isNegativeSemidefinite(position)) throw "out\n";
//            std::cout << objectiveFunction << "\n";
//        if (first) {
        lmi.evaluate(position, settings.B0);
        lmi.evaluateWithoutA0(objectiveFunction, settings.B2);
//        }

        lmi.evaluateWithoutA0(direction, settings.B1);
        MT B2temp = settings.B2 / (-2 * temp);

        // create pencil matrix
        MT AA(2 * matrixDim, 2 * matrixDim);
        MT BB(2 * matrixDim, 2 * matrixDim);

        BB.block(matrixDim, matrixDim, matrixDim, matrixDim) = -1 * settings.B0;
        BB.block(0, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
        BB.block(matrixDim, 0, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
        BB.block(0, 0, matrixDim, matrixDim) = B2temp;

        AA.block(0, matrixDim, matrixDim, matrixDim) = settings.B0;
        AA.block(0, 0, matrixDim, matrixDim) = settings.B1;
        AA.block(matrixDim, 0, matrixDim, matrixDim) = settings.B0;
        AA.block(matrixDim, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);


//        double frobeniusNorm = 0;
//        for (int i=0 ; i<2*matrixDim ; i++)
//            for (int j=0 ; j<2*matrixDim ; j++)
//                frobeniusNorm += BB(i,j) * BB(i,j);
//        frobeniusNorm = std::sqrt(frobeniusNorm);

//        SpMat A =AA.sparseView();
//        SpMat B=BB.sparseView();
//    std::cout << AA << "\n\n\n";
//        std::cout << BB << "\n";

        // Construct matrix operation object using the wrapper classes
//        Spectra::SparseSymMatProd<double> op(A);
        Spectra::DenseSymMatProd<double> op(AA);
//        Spectra::SparseRegularInverse<double> Bop(B);

        Spectra::DenseCholesky<double> Bop(-BB);

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

        for (int i = 0; i < alphas.rows(); i++) {

            if (betas(i) == 0 || alphas(i).imag() != 0)  //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);
//            std::cout << lambda << " e\n";
            if (lambda > 0 && lambda < lambdaMinPositive) {
                lambdaMinPositive = lambda;
                index = i;
            }
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative = lambda;
        }
//        std::cout << lambdaMinPositive << "eval\n";

        Eigen::GeneralizedEigenSolver<MT>::EigenvectorsType eivecs = ges.eigenvectors();
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivec = eivecs.col(index);

        settings.genEigenvector = Point(matrixDim);

        for (int i = 0; i < matrixDim; i++)
            settings.genEigenvector.set_coord(i, eivec(matrixDim + i).real());

//        std::cout << eivecs.col(index)<<"evec\n";
        return lambdaMinPositive;
    }

    template<class Point>
    void compute_reflection(BoundaryOracleBoltzmannHMCSettings<Point> &settings, Point &direction) {
        std::vector<MT> matrices = lmi.getMatrices();
        int dim = matrices.size();
        settings.gradient = Point(dim);

        for (int i = 0; i < dim; i++) {
            settings.gradient.set_coord(i, settings.genEigenvector.dot(
                    (settings.genEigenvector.matrix_left_product(matrices[i]))));
        }

        settings.gradient.normalize();
        settings.gradient = -1 * settings.gradient;

        Point t = ((-2.0 * direction.dot(settings.gradient)) * settings.gradient);
        direction = t + direction;
        settings.gradient = -1 * settings.gradient;
    }

    void compute_reflection(const VT &genEivector, VT &direction, MT &C) {
        VT gradient;
        std::vector<MT> matrices = lmi.getMatrices();
        int dim = matrices.size();
        gradient = VT::Zero(dim);

        for (int i = 0; i < dim; i++) {
            gradient(i) = genEivector.dot((matrices[i]) * genEivector);
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

    template<class Point>
    void compute_reflection(const Point &genEivector, Point &direction) {

        auto t1 = std::chrono::steady_clock::now();

        const std::vector<MT> &matrices = lmi.getMatrices();
        int dim = matrices.size();
        Point gradient(dim);

        for (int i = 0; i < dim; i++) {
            gradient.set_coord(i, genEivector.dot(genEivector.matrix_left_product(matrices[i])));
        }

        gradient = -1 * gradient;
        gradient.normalize();

        gradient = ((-2.0 * direction.dot(gradient)) * gradient);
        direction = direction + gradient;


    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
