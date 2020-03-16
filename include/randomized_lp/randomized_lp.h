#ifndef RANDOMIZED_LP
#define RANDOMIZED_LP

#include <iterator>
#include <vector>
#include <utility>
#include <math.h>
#include <chrono>
#include <unordered_map>
#include "cartesian_geom/cartesian_kernel.h"
#include "vars.h"
#include "hpolytope.h"
#include "interior_point.h"
#include "samplers.h" 

template <typename Parameters, typename Point>
class Randomized_LP{
public:
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

private:
    VT obj;
    MT A; //matrix A
    VT b; // vector b, s.t.: Ax<=b
    bool solved;
    NT objective_value;
    std::vector<NT> variable_values;
    std::unordered_map<unsigned int, NT> lower_bounds;
    std::unordered_map<unsigned int, NT> upper_bounds;

public:
    Randomized_LP(int num_constraints, int num_vars) : A(num_constraints, num_vars), b(num_constraints){
        solved = false;
    }

    // Initialise in the form Ax<=b
    Randomized_LP(const MT &_A, const VT &_b){
        A = _A;
        b = _b;
        solved = false;
    }

    // add constraint as <l, x> <=c 
    bool add_constraint(const VT &_l, int c){
        if(_l.size() != A.cols())
            return false;

        solved = false;

        A.conservativeResize(A.rows()+1, Eigen::NoChange);
        A.row(A.rows()-1) = _l.array();

        b.conservativeResize(A.rows()+1);
        b(A.rows()-1) = c;
    }

    // print the lp model created till now
    bool print(){
        std::cout<<"Objective Function: \n";
        for (size_t i = 0; i < obj.size(); i++){
            if(i==0)
                std::cout<<obj(i)<<"x"<<i;
            else
                std::cout<<" + "<<obj(i)<<"x"<<i;
        }
        std::cout<<"\n";

        std::cout<<"Constraints: \n";
        for (size_t i = 0; i < A.rows(); i++){
            for (size_t j = 0; j<A.cols(); j++){
                if(j==0)
                    std::cout<<A(i,j)<<"x"<<j;
                else
                    std::cout<<" + "<<A(i,j)<<"x"<<j;
            }
            std::cout<<" <= "<<b(i);
            std::cout<<"\n";
        }
        std::cout<<"\n";

        std::cout<<"Bounds: \n";
        for (typename std::unordered_map<unsigned int, NT>::iterator itr = lower_bounds.begin(); itr != lower_bounds.end(); itr++) {
            std::cout<<itr->second<<"<="<<"x"<<itr->first<<"<="<<upper_bounds[itr->first];
        }
        std::cout<<"\n";
    }

    // set column for A
    bool set_column(const VT &c, unsigned int index){
        if(c.rows()!=A.rows())
            return false;

        solved = false;

        A.col(index) = c.array();
        return true;
    }

    // set row for A
    bool set_row(const VT& r, unsigned int index){
        if(r.rows()!=A.cols())
            return false;

        solved = false;

        A.row(index) = r.array();
        return true;
    }

    // set objective function
    bool set_objective(const VT &_obj){
        if(_obj.size()!=A.cols())
            return false;

        solved = false;

        obj = _obj;
        return true;
    }

    // set complete b vector
    bool set_vector(const VT &_b){
        if(b.size() != _b.size())
            return false;

        b = _b;
        return true;
    }

    // set ith index in b vector
    bool set_vector(NT val, unsigned int index){
        b(index) = val;
        return true;
    }

    // set bounds for 'index' numbered variable
    bool set_bounds(unsigned int index, NT lower, NT upper){
        if(index >= A.cols())
            return false;

        lower_bounds[index] = lower;
        upper_bounds[index] = upper;
        return true;
    }

    // get objective function value
    NT get_objective_value(){
        return objective_value;
    }

    // get value of variables after it's solved
    std::vector<NT> get_variables(){
        return variable_values;
    }

    // solve the model created
    bool solve(const Parameters &var, const std::string method="RCP"){
        if(solved == true)
            return true;

        bool res = false;
        if(method == "RCP")
            res = solve_by_cutting_plane(var);
        // else if(method == "SIM_ANN"){
        //     solve_simulated_annealing(var);
        // }
        else
            return false;

        solved = true;
        return res;
    } 

    bool solve_by_cutting_plane(const Parameters &var){
        MT A_curr = A;
        VT b_curr = b;

        unsigned int count = A.rows();
    
        A_curr.conservativeResize(A_curr.rows() + lower_bounds.size() + upper_bounds.size(), Eigen::NoChange);
        b_curr.conservativeResize(b_curr.size() + lower_bounds.size() + upper_bounds.size());

        // Add lower bounds as constraints         
        for (typename std::unordered_map<unsigned int, NT>::iterator itr = lower_bounds.begin(); itr != lower_bounds.end(); itr++, count++) { 
            unsigned int var_num = itr->first;
            NT bnd_val = itr->second;
            VT temp(A_curr.cols());
            for(int i=0; i<temp.size(); i++)
                temp(i)=0;
            temp(var_num) = -1;
            A_curr.row(count) = temp.array();
            b_curr(count) = -bnd_val;
        }
        
        // Add upper bounds as constraints
        for (typename std::unordered_map<unsigned int, NT>::iterator itr = upper_bounds.begin(); itr != upper_bounds.end(); itr++, count++) { 
            unsigned int var_num = itr->first;
            NT bnd_val = itr->second;
            VT temp(A_curr.cols());
            for(int i=0; i<temp.size(); i++)
                temp(i)=0;
            temp(var_num) = 1;
            A_curr.row(count) = temp.array();
            b_curr(count) = bnd_val;
        }

        // Using objective fxn as a point for dot product
        std::vector<NT> tmp1(obj.data(), obj.data() + obj.size());
        Point obj_c_pt(obj.rows(), tmp1);

        unsigned int max_iterations = 100;
        double threshold = 0.001;
        double N_k_init = 2.2*log(1/threshold) + 1.1 + 5;
        NT min_val = std::numeric_limits<NT>::max();
        Point p;
        HPolytope<Point> pol;
        
        A_curr.conservativeResize(A_curr.rows()+1, Eigen::NoChange);
        b_curr.conservativeResize(b_curr.rows()+1, Eigen::NoChange);
        count = A_curr.rows()-1;
        A_curr.row(count) = obj.array();
        b_curr(count) = 0;
        pol.init(A_curr.cols(), A_curr, b_curr);
        
        for (unsigned int k = 0; k < max_iterations; k++){
            // unsigned int N_k = (N_k_init + 0.505*log(k+1));
            unsigned int N_k = 16;

            pol.set_vec(b_curr);
            // pol.print();
            p = getInteriorPoint<Point, NT>(pol);
            if(pol.is_in(p)==0)
                break;
            // p.print();
            // std::cout<<"Var n: "<<var.n<<", k: "<<k<<std::endl;

            std::vector<Point> randPoints; //ds for storing rand points
            rand_point_generator(pol, p, N_k, 1000, randPoints, var);
            for (typename std::vector<Point>::iterator it = randPoints.begin(); it != randPoints.end(); ++it){
                // evaluating <obj,pt>
                NT val = obj_c_pt.dot(*it);
                // (*it).print();
                // std::cout<<it->get_coeffs().back()<<"--"<<val<<"\n";
                if(val < min_val){
                    min_val = val;
                }
            }

            b_curr(count) = min_val;
        }

        variable_values = p.get_coeffs();
        objective_value = min_val;

        return true;
    }

    // bool solve_by_simulated_annealing(const Parameters &var){

    // }
};

#endif