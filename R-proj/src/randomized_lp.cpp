#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "randomized_lp.h"

//' Solving a LP problem using Randomized algorithms.
//'
//' Currently, only 2 algorithms are implemented - Randomized Cutting plane and Simulated Annealing.
//'
//' @param P. A convex H Polytope, it is the feasible region of the LP problem (\eqn{Ax<=b}).
//' @param obj. A vector for the coefficients of the objective function (\eqn{min c^{T} x}).
//' @param bounds optional. A list that contains the bound of the variables (the default is complete Real space), as follows:
//' \itemize{
//' \item{\code{indices}}{A vector containing the variable indices(0 indexed) for which the bounds have to be set.}
//' \item{\code{lower}}{A vector containing the value of lower bounds for all the variables specified in indices.}
//' \item{\code{upper}}{A vector containing the value of upper bounds for all the variables specified in indices.}
//' }
//' @param algo Optional. An unsigned integer that declares which algorithm, as follows:
//' \itemize{
//' \item{\code{0} }{ Use the Randomized Cutting Plane algorithm (RCP).}
//' \item{\code{1} }{ Use the Simulated Annealing algorithm (SIM_ANN).}
//' }
//' @param verbose Optional. A boolean parameter for printing out the LP program formed.
//'
//' @references \cite{ Dabbene, Fabrizio, Pavel S. Shcherbakov, and Boris T. Polyak.,
//' \dQuote{ A randomized cutting plane method with probabilistic geometric convergence,} \emph{SIAM Journal on Optimization 20.6,} (2010): 3185-3207.},
//' @references \cite{Adam Tauman Kalai, Santosh Vempala,
//' \dQuote{Simulated Annealing for Convex Optimization,} \emph{Mathematics of Operations Research Vol. 31, No. 2,} 2006.}
//'
//'
//' @return A list containing the value of the objective function and value of all variables.
//' @examples
//' # computing Chebychev ball for a H-polytope (3d cube)
//' P <- gen_cube(3, 'H')
//' row_norm <- sqrt(rowSums((P$A)^2))
//' P$A <- cbind(P$A, row_norm)
//' var_bounds <- list("indices"=c(3), "lower"=c(0), "upper"=c(1000))
//' randomized_lp_solver(P, obj=c(0,0,0,-1), bounds=var_bounds, algo=1, verbose=TRUE)
//' @export
// [[Rcpp::export]]
Rcpp::List randomized_lp_solver(Rcpp::Reference P, Rcpp::NumericVector obj, Rcpp::Nullable<Rcpp::List> bounds = R_NilValue, unsigned int algo = 0, bool verbose=false){
    typedef double NT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    // Setup the parameters
    unsigned int n = P.field("dimension");
    unsigned int walk_len=10 + n/10;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001;
    int rnum = 100; // will be configured by the algorithm
    RNGType rng(std::time(0));
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars<NT,RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,false, true,false);

    MT A = Rcpp::as<MT>(P.field("A"));
    VT b = Rcpp::as<VT>(P.field("b"));
    VT obj_fxn = Rcpp::as<VT>(obj);

    if(obj_fxn.size() != A.cols())
        throw Rcpp::exception("Invalid objective vector, the size of objective vector does not match with Polytope's dimension");

    std::string method = "";
    if(algo >= 2)
        throw Rcpp::exception("Invalid algo type, currently volesti supports only 2 algorithms for randomized lp solver");
    else if(algo == 0)
        method = "RCP";
    else
        method = "SIM_ANN";

    Randomized_LP<vars<NT,RNGType>, Point> lp_solver(A, b);
    lp_solver.set_objective(obj_fxn);

    if(bounds.isNotNull()){
        if((!Rcpp::as<Rcpp::List>(bounds).containsElementNamed("lower")) || (!Rcpp::as<Rcpp::List>(bounds).containsElementNamed("indices")) || (!Rcpp::as<Rcpp::List>(bounds).containsElementNamed("upper"))){
            throw Rcpp::exception("Invalid bounds argument");
        }

        Rcpp::NumericVector bounds_indices = Rcpp::as<Rcpp::NumericVector>(Rcpp::as<Rcpp::List>(bounds)["indices"]);
        Rcpp::NumericVector lower_bounds = Rcpp::as<Rcpp::NumericVector>(Rcpp::as<Rcpp::List>(bounds)["lower"]);
        Rcpp::NumericVector upper_bounds = Rcpp::as<Rcpp::NumericVector>(Rcpp::as<Rcpp::List>(bounds)["upper"]);

        if(bounds_indices.length() != lower_bounds.length() || bounds_indices.length() != upper_bounds.length()){
            throw Rcpp::exception("Invalid bounds argument, length of all 3 vectors should be same");
        }

        for(int i=0; i<bounds_indices.length(); i++) {
            if(!lp_solver.set_bounds(bounds_indices[i], lower_bounds[i], upper_bounds[i]))
                throw Rcpp::exception("Invalid bound!!");
        }
    }

    if(verbose){
        lp_solver.print();
        Rcpp::Rcout<<"\n";
    }

    if(!lp_solver.solve(var, method))
        throw Rcpp::exception("The solver could not solve the provided LP problem");

    std::vector<NT> variable_values = lp_solver.get_variables();
    NT obj_value = lp_solver.get_objective_value();

    Rcpp::List res;
    res["objective_value"] = obj_value;
    res["variables"] = Rcpp::wrap(variable_values);
    return res;
}
