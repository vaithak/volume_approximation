
#ifndef ZONO_OVER_APPROX_H
#define ZONO_OVER_APPROX_H

template<class Hpolytope, class Zonotope>
void zono_approx (Zonotope &Z, bool hp) {

    typedef typename Zonotope::NT NT;
    typedef typename Zonotope::PolytopePoint Point;
    typedef typename Zonotope::MT MT;
    typedef typename Zonotope::VT VT;
    typedef boost::mt19937 RNGType;

    int n = Z.dimension(), k = Z.num_of_generators();
    NT e = 0.1, delta = -1.0, lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0, alpha = 0.2, diam = -1.0;
    int win_len = 150, NN = 125, nu =10, walkL = 1;
    bool ball_walk = false, verbose = false, cdhr = false, rdhr = false, billiard = true, round = false, win2 = false, hpoly = false;


    NT ratio = 0.0, vol_red = 0.0, vol = 0.0;

    MT X(n, 2*k);
    X << Z.get_mat().transpose(), -Z.get_mat().transpose();
    Eigen::JacobiSVD<MT> svd(X*X.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    MT G(k, 2*n);
    G << Z.get_mat()*svd.matrixU(), Z.get_mat()*svd.matrixU();
    VT Gred_ii = G.transpose().cwiseAbs().rowwise().sum();
    MT A(n, 2*n);
    A << -MT::Identity(n,n), MT::Identity(n,n);
    MT Mat(2*n, n+1);

    Mat << Gred_ii, A.transpose()*svd.matrixU().transpose();

    Hpolytope HP;
    HP.init(n, A.transpose()*svd.matrixU().transpose(), Gred_ii);

    //if (true) {
        vol_red = std::abs(svd.matrixU().determinant());
        for (int i = 0; i < n; ++i) {
            vol_red *= 2.0 * Gred_ii(i);
        }

        //e = error;
        hpoly = hp;

        if (billiard && diam < 0.0) Z.comp_diam(diam);

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        // the random engine with this seed
        typedef boost::mt19937    RNGType;
        RNGType rng(seed);
        boost::random::uniform_real_distribution<>(urdist);
        boost::random::uniform_real_distribution<> urdist1(-1,1);

        std::pair<Point,NT> InnerB;
        InnerB = Z.ComputeInnerBall();

        vars<NT, RNGType> var(1, n, walkL, 1, 0.0, e, 0, 0.0, 0, InnerB.second, diam, rng,
                               urdist, urdist1, delta, false, false, round, false, false, ball_walk, cdhr,rdhr, billiard,
                               0.0, 0.0, 0.0);
        vars_ban <NT> var_ban(lb, ub, p, 0.0, alpha, win_len, NN, nu, win2);


        NT nballs;
        double tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        if (!hpoly) {
            vol = volesti_ball_ann(Z, var, var_ban, InnerB, nballs);
        } else {
            vars_g<NT, RNGType> varg(n, 1, 1000 + NT(n * n) / 2.0, 6*n*n+500, 1, e, InnerB.second, rng, 2.0, 0.1,
                                     1.0 - 1.0 / (NT(n)), delta, false, false, false, false, false, false,
                                     false, true, false, 0.0, 0.0);
            vol = vol_hzono<Hpolytope> (Z, var, var_ban, varg, InnerB, nballs);
        }
        double tstop = (double)clock()/(double)CLOCKS_PER_SEC;
        ratio = std::pow(vol_red / vol, 1.0/NT(n));
    //}

    std::cout<<"estimated volume = "<< vol<<"\nnumber of balls in MMC = "<<nballs<<"\nnumber of Boundary Oracle Calls = "<<var.BoundCalls<<"\ntime = "<<tstop-tstart<<" sec"<<std::endl;
    std::cout<<"volume of PCA overapproximation polytope = "<<vol_red<<"\nratio of fitness = "<<ratio<<std::endl;


}

#endif // ZONO_OVER_APPROX_H
