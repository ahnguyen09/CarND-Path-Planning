/*
Implementing JMT in a single header file like spline.h 
*/

#ifndef JMT_H
#define JMT_H

#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include <vector>
#include "Dense"

namespace jmt
{

std::vector<double> JMT(std::vector< double> start, std::vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT 
    an array of length 6, each value corresponding to a coefficent in the polynomial 
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */
    
    double T_2 = T*T;
	double T_3 = T_2*T;
	double T_4 = T_3*T;
	double T_5 = T_4*T;
	
	Eigen::MatrixXd three_coeff = Eigen::MatrixXd(3, 3);

    three_coeff << T_3,    T_4,   T_5,
                   3*T_2,    4*T_3,   5*T_4,
                   6*T,    12*T_2,   20*T_3;
    
    Eigen::VectorXd three_conditions = Eigen::VectorXd(3);
    
    three_conditions << end[0] - (start[0] + start[1]*T + 0.5*start[2]*T_2),
                        end[1] - (start[1] + start[2]*T),
                        end[2] - start[2];

	Eigen::VectorXd ans = three_coeff.colPivHouseholderQr().solve(three_conditions);


    return {start[0],start[1],0.5*start[2],ans(0),ans(1),ans(2)}; 
}

} //jmt namespace

#endif /* JMT_H */