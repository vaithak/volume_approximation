#ifndef SOLVE_RANDOMIZED_LP_H
#define SOLVE_EANDOMIZED_LP_H

#include <stdio.h>
#include <cmath>
#include <exception>
#include "samplers.h"
#include "randomized_lp.h"


// compute the chebychev ball of an H-polytope described by a dxd matrix A and  d-dimensional vector b, s.t.: Ax<=b
template <typename NT, typename Point, typename MT, typename VT, class Parameters>
std::pair<Point,NT> LPChebychevBall(const MT &A, const VT &b, Parameters &var){
    int d = A.cols();
    int Ncol=d+1, m=A.rows();
    int *colno = NULL;
    Randomized_LP<Parameters, Point> lp(m, Ncol);
    std::pair<Point,NT> exception_pair(Point(1),-1.0);
    //std::cout<<"Hi, reached a"<<std::endl;
    try{
        if(!lp.set_vector(b)) throw false;
    }
    catch (bool e){
        #ifdef VOLESTI_DEBUG
        std::cout<<"Could not add vector b for the Linear Program for chebychev center "<<e<<std::endl;
        #endif
        return exception_pair;
    }
    //std::cout<<"Hi, reached 2"<<std::endl;
    for (int i = 0; i < m; ++i) {
        /* construct all rows */
        VT ith_row(Ncol);
        
        for(int j=0; j<d; j++){
            ith_row(j) = A(i,j);
        }
        ith_row(d) = A.row(i).norm();
        
        /* add the row to lpsolve */
        try {
            if(!lp.set_row(ith_row, i)) throw false;
        }
        catch (bool e)
        {
            #ifdef VOLESTI_DEBUG
            std::cout<<"Could not define constriants for the Linear Program for chebychev center "<<e<<std::endl;
            #endif
            return exception_pair;
        }
    }
    //std::cout<<"Hi, reached 3"<<std::endl;

    NT infinite = std::numeric_limits<NT>::max();
    lp.set_bounds(d, 0.0, 1000);

	// set the objective function
    VT obj(Ncol);
    for(int i=0; i<d; i++){
        obj(i)=0;
    }
    //std::cout<<"Hi, reached 4"<<std::endl;
    // because we have to minimize -r
    obj(d) = -1;
    try{
        if (!lp.set_objective(obj)) throw false;
    }
    catch (bool e)
    {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Could not define objective function for the Linear Program for chebychev center "<<e<<std::endl;
        #endif
        return exception_pair;
    }

    //std::cout<<"Hi, reached 5"<<std::endl;
    // lp.print();
    var.n += 1;
    /* Now let lpsolve calculate a solution */
    try
    {
        if (!lp.solve(var)) throw false;
    }
    catch (bool e)
    {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Could not solve the Linear Program for chebychev center "<<e<<std::endl;
        #endif
        return exception_pair;
    }

    //std::cout<<"Hi, reached 6"<<std::endl;
    std::pair<Point,NT> res;
    std::vector<NT> temp_p = lp.get_variables();
    Point xc( d , temp_p.begin() , --temp_p.end());
    NT r= -lp.get_objective_value();
    //std::cout<<"Hi, reached 7"<<std::endl;
    res = std::pair<Point,NT> (xc,r);

    return res;
}

#endif