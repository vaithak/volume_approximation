//
// Created by panagiotis on 9/7/2019.
//

#ifndef VOLESTI_SDP_GENERATOR_H
#define VOLESTI_SDP_GENERATOR_H

#include "spectrahedron.h"
#include "Eigen"
#include "sdp_problem.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
typedef boost::mt19937 RNGType;


void randomMatrixGOE(MT& M) {
    unsigned m = M.rows();
    boost::normal_distribution<> rdist(0,1);
    unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();//4 if fixed for debug
    RNGType rng(seed);

    for (unsigned int i=0; i<m; i++) {
        for (unsigned int j=0; j<m; j++) {
            M(i,j) = rdist(rng);
        }
    }
}
// m must be even

template <class Point>
optimization::sdp_problem<Point> generateSDP(int n, int m) {
    VT obj(n);
    srand (time(NULL));
    
    if (m % 2 != 0) throw 1;

    for (int i = 0; i < n; i++) {
        obj(i) = -5 + ((double)rand() / RAND_MAX)*30;
    }

    MT ones(m, m);
    for (int i=0 ; i<ones.rows() ; i++)
        for (int j=0 ; j<ones.cols() ; j++)
            ones(i, j) = 1;

    MT M = 2* Eigen::MatrixXd::Random(m,m) - ones;
//    M.resize(m,m);
//    randomMatrixGOE(M);



    MT I = Eigen::MatrixXd::Identity(m, m);
    std::vector<MT> matrices(n + 1);
    matrices[0] = -(M * M.transpose()) - I;

//    matrices[0] = I;//(M + M.transpose())/2;

    ones.resize(m/2, m/2);
    for (int i=0 ; i<ones.rows() ; i++)
        for (int j=0 ; j<ones.cols() ; j++)
            ones(i, j) = 1;


    for (int i=1 ; i<=n ; i++) {
        M =  2* Eigen::MatrixXd::Random(m/2,m/2) - ones;
        MT MM = M.selfadjointView<Eigen::Upper>();

        for (int at = 0 ; at<m/2 ; at++)
            MM(at,at) -= M(at,at);

        MT A;
        A.setZero(m, m);
        A.topLeftCorner(m/2, m/2) = MM;
        A.bottomRightCorner(m/2, m/2) = -MM;
        matrices[i] = A;

//        M.resize(m,m);
//        randomMatrixGOE(M);
//        matrices[i] = (M + M.transpose())/2- ones;
    }

    LMI lmi(matrices);
    Spectrahedron spectrahedron(lmi);

    return optimization::sdp_problem<Point>(spectrahedron, obj);
}

#endif //VOLESTI_SDP_GENERATOR_H
