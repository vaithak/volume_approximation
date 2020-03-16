#ifndef INTERIOR_POINT_H
#define INTERIOR_POINT_H

#include "hpolytope.h"

template<typename NT, typename MT, typename VT>
NT fxn(NT s, const VT& x, const MT& A, const VT& b, double m) {
    NT ret = s;

    for (int i=0 ; i<A.rows(); i++) {
        NT temp = s - A.row(i).dot(x) + b(i);
        if (temp <= 0.00001)
            return std::numeric_limits<NT>::max();
        ret -=  m*log(temp);
    }
    return ret;
}

template<typename NT, typename VT, typename MT>
double backtrackingLineSearch(NT s, const VT &grad, const NT &s_grad, const MT &A, const VT &b, VT &x, double m) {
    double step_length= 1;
    double alpha = 0.003;
    double beta = 0.7;

    VT eval = x - step_length * grad;
    NT eval_s = s - step_length*s_grad;

    // Armijo-Goldstein condition
    while (fxn(eval_s, eval, A, b, m) >= (fxn(s, x, A, b, m) - alpha*step_length * (s_grad*s_grad + grad.dot(grad)))) {
        step_length *= beta;
        eval = x - step_length * grad;
        eval_s = s - step_length*s_grad;

        if(step_length < 1e-19)
            return step_length;
    }

    return step_length;
}

template<typename NT, typename MT, typename VT>
void gradientDescent(const MT &A, const VT &b, NT &s, VT &x, const double m) {
    long dim = A.cols();
    VT denominators(A.rows());  // the denominators from the derivatives of log are common at each iteration
    VT grad(dim);
    NT s_grad;

    int max_iter = 1000;
    do {
        // compute the gradient
        for (int i = 0; i < denominators.size(); i++) {
            denominators(i) =  m / (s - A.row(i).dot(x) + b(i));
        }

        s_grad = 1;
        for (int i = 0; i < denominators.size(); i++)
            s_grad -=  denominators(i);

        for (int j = 0; j < dim; j++) {
            grad(j)=0;
            for (int i = 0; i < A.rows(); i++)
                grad(j) += A(i, j) * denominators(i);
        }

        // we can also exit if s < 0. It suffices for our problem
        if (sqrt(s_grad * s_grad + grad.dot(grad)) < 0.001 || s < 0) 
            break;

        // compute step length
        NT step_length = backtrackingLineSearch(s, grad, s_grad, A, b, x, m);
        // std::cout<<grad.norm()<<", S grad: "<<s_grad<<", Norm grad:  "<<sqrt(s_grad * s_grad + grad.dot(grad))<<", step: "<<step_length<<", s: "<<s<<", m: "<<m<<std::endl;
        // update point
        x -= step_length * grad;
        s -= step_length * s_grad;

    } while (max_iter--);

    if(max_iter <= 0){
        #ifdef VOLESTI_DEBUG
        // std::cout<<"Max iterations exceeded in gradient descent\n";
        #endif
    }
}


/**
 * Returns a feasible point in the polytope. It uses barrier method to solve the phase I problem to get a feasible point
 *
 * minimize s
 * s.t. f_i(x) <= s
 *
 * The resulting problem is:
 *
 * s - m(log(s - f_1) + log(s - f_2) + ... ), m going to 0
 * If failed at finding a point, throw exception
 *
 * @tparam Point
 * @tparam NT
 * @param polytope
 * @return an feasible point
 */
template<class Point, typename NT>
Point getInteriorPoint(HPolytope<Point>& polytope) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;


    MT A = polytope.get_mat();
    VT b = polytope.get_vec();
    unsigned int dim = polytope.dimension();

    NT s;
    VT x(dim);

    // the initial value for x will be 1 and for s = max{f_i(x)} + 1
    for (int j=0 ; j<dim ; j++)
        x(j) = 1;

    s = A.row(0).dot(x) - b(0);

    for (int i=1 ; i<A.rows() ; i++) {
        NT _s = A.row(i).dot(x) - b(i);
        if (_s > s)
            s = _s;
    }
    s += 1;

    unsigned int steps = 10;
    double m = 1;

    while (steps > 0) {
        gradientDescent(A, b, s, x, m);
        // std::cout<<"Break:-> "<<"m: "<<m<<"s: "<<s<<std::endl;
        if (s < 0) break; // it suffices to get s below 0
        steps--;
        m /= 10;
    }

    // if (s>0){
        // throw std::string("No internal point was found in the polytope");
    // }

    std::vector<NT> tmp1(x.data(), x.data() + x.size());
    Point x_pt(dim, tmp1);
    return x_pt;
}

#endif //INTERIOR_POINT_H