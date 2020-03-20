#ifndef HELPER_H
#define HELPER_H

#include "hpolytope.h"
#include "samplers.h"

template <typename NT, typename MT, typename Point>
static void calculate_covariance_matrix(std::vector<Point> &points, MT &cov_mat){
    if(points.size() == 0)
        return;
        
    MT temp_matrix(points.size(), points[0].dimension());

    for(unsigned int i = 0; i<temp_matrix.rows(); ++i){
        std::vector<NT> temp_vec = points[i].get_coeffs();
        temp_matrix.row(i) = Eigen::Map<Eigen::VectorXd>(points[i].get_coeffs().data() ,points[i].dimension());
    }

    temp_matrix = temp_matrix.rowwise() - temp_matrix.colwise().mean();
    cov_mat = temp_matrix.transpose() * temp_matrix;

    if(points[0].dimension() != 1)
        cov_mat = cov_mat/(double)(points[0].dimension() - 1);
}

// radius of ball containing the polytope centered at x
template <typename NT, typename Point>
std::pair<NT, NT> calculate_max_min_dist(HPolytope<Point> pol , Point x){
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    NT res_min = std::numeric_limits<NT>::max();
    NT res_max = std::numeric_limits<NT>::min();
    unsigned int num_planes = pol.num_of_hyperplanes();
    Eigen::Map<Eigen::VectorXd> x_pt(x.get_coeffs().data(), x.dimension());
    MT A = pol.get_mat();
    VT b = pol.get_vec();

    for (int i=0; i<num_planes; i++){
        if(A.row(i).norm()==0)  continue;
        
        // distance of x from ith hyperplane
        NT temp = abs(A.row(i).dot(x_pt) - b(i))/A.row(i).norm();

        if(temp > res_max)  res_max = temp;
        if(temp < res_min)  res_min = temp;
    }

    return {res_min, res_max};
}

// Density function for exponential distribution
double dexp(double x, double lambda){
    return lambda*exp(-lambda*x);
}

// Cumulative distribution function for exponential distribution
double pexp(double x, double lambda){
    return 1 - exp(- lambda*x);
}

// Quantile function for exponential distribution
double qexp(double p, double lambda){
    return -log(1 - p) / lambda;
}

// Produce a randomly generated number from a truncated exponential distribution using inverse CDF method
template <class RNGType>
double texp(double lambda, double a, double b, RNGType& rng) {
    boost::random::uniform_real_distribution<> urdist(0, 1);
    double u = urdist(rng);
    double cdfA = pexp(a, lambda);
    double cdfB = pexp(b, lambda);

    if((cdfA-cdfB)< 1e-10)
        return a;

    return qexp(cdfA + u*(cdfB - cdfA), lambda);
}

// generating a point for simulated annealing with hit and run random walk with boltzmann function
template <class Polytope, class Point, class Parameters, typename MT>
Point boltzmann_hit_and_run_annealing(Point &start_pt, Polytope &P, Point& cost_fxn, double Temperature, const MT& cholesky_decomp, const Parameters &var){
    typedef typename Parameters::RNGType RNGType;
    RNGType &rng = var.rng;
    unsigned int dim = start_pt.dimension();
    typedef typename Point::FT NT;

    // generate direction with distrubtion N(0, I), I - identity matrix of (dim x dim)
    Point temp_dir = get_direction<RNGType, Point, NT>(dim);
    std::vector<NT> temp_dir_vec = temp_dir.get_coeffs();

    // generate direction with distrubtion N(0, V), V - cholesky*choleskyT of (dim x dim)
    Eigen::VectorXd dir_vec = cholesky_decomp*Eigen::Map<Eigen::VectorXd>(temp_dir.get_coeffs().data(), temp_dir.dimension());
    Point dir(dim, std::vector<NT>(dir_vec.data(), dir_vec.data() + dir_vec.size()));
    
    std::pair<NT, NT> intersect_lambdas = P.line_intersect(start_pt, dir);
    Point pol_intersect_pt_1 = (intersect_lambdas.first * dir) + start_pt;
    Point pol_intersect_pt_2 = (intersect_lambdas.second * dir) + start_pt;
    double cost_1 = cost_fxn.dot(pol_intersect_pt_1);
    double cost_2 = cost_fxn.dot(pol_intersect_pt_2);

    NT lambda = texp(abs(cost_1 - cost_2)/Temperature, 0.001, 0.999, rng);
    Point ext_dir = (pol_intersect_pt_1 - pol_intersect_pt_2);
    
    if(cost_1 > cost_2)
        return (lambda * ext_dir)  + pol_intersect_pt_2;

    return (-lambda * ext_dir) + pol_intersect_pt_1;
}

#endif