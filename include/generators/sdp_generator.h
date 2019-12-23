//
// Created by panagiotis on 9/7/2019.
//

#ifndef VOLESTI_SDP_GENERATOR_H
#define VOLESTI_SDP_GENERATOR_H

//#include "spectrahedron.h"
//#include "Eigen"
//#include "sdp_problem.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
typedef boost::mt19937 RNGType;


template <class MT>
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

template <class LMI, class Spectrahedron, class Point>
Spectrahedron generateSDP(int n, int m) {

    typedef typename Spectrahedron::VT VT;
    typedef typename Spectrahedron::MT MT;
    typedef typename Point::FT NT;

    MT ones = MT::Ones(m, m);
    MT M = 2* Eigen::MatrixXd::Random(m,m) - ones;

    MT I = Eigen::MatrixXd::Identity(m, m);
    std::vector<MT> matrices(n + 1);
    matrices[0] = -(M * M.transpose()) - I;


    MT ones2 = MT::Ones(m/2, m/2);
    MT MM(m/2, m/2);

    for (int i=1 ; i<=n ; i++) {
        MM =  2 * MT::Random(m/2, m/2) - ones2;
        MM = MM + MM.transpose();
        MT A;
        A.setZero(m, m);

        for (int j = 0; j < m/2; ++j) {
            for (int k = 0; k < m/2; ++k) {
                A(j,k) = MM(j,k);
                A(j+m/2, k+m/2) = -MM(j+m/2, k+m/2);
            }
        }
        std::cout<<"A = "<<A<<"\n"<<std::endl;
        matrices[i] = A;
    }

    LMI lmi(matrices);
    Spectrahedron spectrahedron(lmi);
    return spectrahedron;

    //return optimization::sdp_problem<Point>(spectrahedron, obj);
}

#endif //VOLESTI_SDP_GENERATOR_H
