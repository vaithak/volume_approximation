#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "misc.h"
#include "known_polytope_generators.h"
#include "solve_randomized_lp.h"


template <typename NT, class RNGType, class Polytope>
void cheb_test(Polytope &P, NT expected, NT tolerance=0.0001)
{
    typedef typename Polytope::PolytopePoint Point;

    // Setup the parameters
    int n = P.dimension();
    int walk_len=10 + n;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001;
    int rnum = 5*n;
    RNGType rng(std::time(0));
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars<NT,RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,false,true,false);

    //Compute chebychev ball//
    //std::cout << "\n--- Testing Chebchev ball computation of " << f << std::endl;
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    MT A = P.get_mat();
    VT b = P.get_vec();
    std::pair<Point,NT> CheBall = LPChebychevBall<NT, Point>(A, b, var, method="SIM_ANN");
    double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;


    std::cout<<"\nradius is: "<<CheBall.second<<std::endl;
    std::cout << "Chebychev time = " << tstop1 - tstart1 << std::endl;

    NT error = std::abs(CheBall.second-expected)/expected;

    CHECK(error < tolerance);
}

template <typename NT>
void call_test_cube() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing Chebchev ball computation of H-cube10" << std::endl;
    P = gen_cube<Hpolytope>(10, false);
    cheb_test<NT, RNGType>(P, 1.0, 0.05);

    std::cout << "\n--- Testing Chebchev ball computation of H-cube20" << std::endl;
    P = gen_cube<Hpolytope>(20, false);
    cheb_test<NT, RNGType>(P, 1.0, 0.05);
}

template <typename NT>
void call_test_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex10" << std::endl;
    P = gen_simplex<Hpolytope>(10, false);
    cheb_test<NT, RNGType>(P, 0.0759747, 0.05);

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex20" << std::endl;
    P = gen_simplex<Hpolytope>(20, false);
    cheb_test<NT, RNGType>(P, 0.0408628, 0.05);

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex30" << std::endl;
    P = gen_simplex<Hpolytope>(30, false);
    cheb_test<NT, RNGType>(P, 0.0281871, 0.05);

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex40" << std::endl;
    P = gen_simplex<Hpolytope>(40, false);
    cheb_test<NT, RNGType>(P, 0.0215868, 0.05);
}

TEST_CASE("rcp_cheb_cube") {
    call_test_cube<double>();
}

TEST_CASE("rcp_cheb_simplex") {
    call_test_simplex<double>();
}
