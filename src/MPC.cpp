#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
//#include "cppad/cppad.hpp"
//#include "cppad/ipopt/solve.hpp"
#include <iostream>

#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.2;


// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.


// init counters
const int x_start = 0;
const int y_start = x_start + N;
const int psi_start = y_start + N;
const int v_start = psi_start + N;
const int cte_start = v_start + N;
const int epsi_start = cte_start + N;
const int delta_start = epsi_start + N;
const int a_start = delta_start + N - 1;

//
// MPC class definition implementation.
//
MPC::MPC() {
	// variables for green line in sim
	glob_x.resize(N - 1);
	glob_y.resize(N - 1); 
	Lf = 2.67;
}

class FG_eval {
public:
	// Fitted polynomial coefficients
	Eigen::VectorXd coeffs;
	FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

	typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
	void operator()(ADvector& fg, const ADvector& vars) {
		// TODO: implement MPC
		// `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
		// NOTE: You'll probably go back and forth between this function and
		// the Solver function below.
		
		//for (int i = 0; i < N; i++) {
		//	std::cout << "vars[v_start+i]" << vars[v_start+i] << std::endl;
		//}

		// reset
		fg[0] = 0;
		for (int i = 0; i < N; i++) {
			// add cross track error [4]
			//cout << "vars[cte_start + i]   -  " << vars[cte_start + i] << endl;
			//cout << "vars[" << i  << "]   -  " << vars[i] << endl;

			fg[0] += CppAD::pow(vars[cte_start + i], 2);
			fg[0] += CppAD::pow(vars[epsi_start+ i], 2);
			fg[0] += CppAD::pow(vars[v_start+ i]-40, 2);
		}

		// minimize actuators
		for (int i = 0; i < N - 1; i++) {
			fg[0] += CppAD::pow(vars[i+delta_start], 2);
			fg[0] += CppAD::pow(vars[i+a_start], 2);
		}

		// minimize gap
		for (int i = 0; i < N - 2; i++) {
			fg[0] += 1000.0 * CppAD::pow( vars[delta_start+i+1] - vars[delta_start+i], 2);
			fg[0] += 1000.0 * CppAD::pow(vars[i+7*N] - vars[i+7*N], 2);
			//fg[0] += CppAD::pow( vars[delta_start+i+1] - vars[delta_start+i], 2);
			//fg[0] += CppAD::pow(vars[i+7*N] - vars[i+7*N], 2);
		}


		fg[1 + x_start]   = vars[x_start];
		fg[1 + y_start]   = vars[y_start];
		fg[1 + psi_start] = vars[psi_start];
		fg[1 + v_start]   = vars[v_start];
		fg[1 + cte_start] = vars[cte_start];
		fg[1 + epsi_start]= vars[epsi_start];

		// model from the class
		for (int i = 0; i < N - 1; i++) {
			fg[2+i + x_start]   = vars[1+i + x_start]   - (vars[i + x_start]   + vars[i + v_start] * CppAD::cos(vars[i + psi_start]) * dt);
			fg[2+i + y_start]   = vars[1+i + y_start]   - (vars[i + y_start]   + vars[i + v_start] * CppAD::sin(vars[i + psi_start]) * dt);
			fg[2+i + psi_start] = vars[1+i + psi_start] - (vars[i + psi_start] + vars[i + v_start] * vars[i + delta_start] / 2.67    * dt);
			fg[2+i + v_start]   = vars[1+i + v_start]   - (vars[i + v_start]   + vars[i + a_start]                                   * dt);


			AD<double> f0 = coeffs[0] + coeffs[1] * vars[i + x_start] + coeffs[2] * pow(vars[i + x_start],2) + coeffs[3] * pow(vars[i + x_start], 3);
			AD<double> side0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * vars[i + x_start] + 3 * coeffs[3] * pow(vars[i + x_start], 2));

			fg[2+i + cte_start]  = vars[1+i + cte_start]  - ( (f0 - vars[i + y_start]     ) + vars[i + v_start] * CppAD::sin(vars[i + epsi_start]) * dt);
			fg[2+i + epsi_start] = vars[1+i + epsi_start] - ( (vars[i + psi_start] - side0) + vars[i + v_start] * vars[i + delta_start] / 2.67     * dt);
		}
		//std::cout << "here 7\n";
	}

	
};


MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
	bool ok = true;
	size_t i;
	typedef CPPAD_TESTVECTOR(double) Dvector;
	
	// TODO: Set the number of model variables (includes both states and inputs).
	// For example: If the state is a 4 element vector, the actuators is a 2
	// element vector and there are 10 timesteps. The number of variables is:
	//
	// 4 * 10 + 2 * 9
	
	// 6 variables
	// N steps
	// 2 actuators
	// Start with 0 (about '-1')
	size_t n_vars = 6*N+2*(N-1);
	// TODO: Set the number of constraints
	size_t n_constraints = 6*N;


	
	// Initial value of the independent variables.
	// SHOULD BE 0 besides initial state.
	Dvector vars(n_vars);
	for (int i = 0; i < n_vars; i++) {
		//reset
		vars[i] = 0;
	}
	// init state
	vars[x_start]     = state[0];
	vars[y_start]     = state[1];
	vars[psi_start]   = state[2];
	vars[v_start]     = state[3];
	vars[cte_start]   = state[4];
	vars[epsi_start]  = state[5];

	//cout << "vars[epsi_start]   - " << vars[epsi_start] << endl;

	Dvector vars_lowerbound(n_vars);
	Dvector vars_upperbound(n_vars);
	// TODO: Set lower and upper limits for variables.

	for (int i = 0; i < delta_start; i++) {
		vars_lowerbound[i] = -1.0e19;
		vars_upperbound[i] = 1.0e19;
	}
	//-25 to 25 degrees
	for (int i = delta_start; i < a_start; i++) {
		vars_lowerbound[i] = -0.436332;
		vars_upperbound[i] = 0.436332;
	}

	for (int i = a_start; i < n_vars; i++) {
		// from -1 to 1
		vars_lowerbound[i] = -1.0;
		vars_upperbound[i] = 1.0;
	}


	// Lower and upper limits for the constraints
	// Should be 0 besides initial state.
	Dvector constraints_lowerbound(n_constraints);
	Dvector constraints_upperbound(n_constraints);
	for (int i = 0; i < n_constraints; i++) {
		constraints_lowerbound[i] = 0;
		constraints_upperbound[i] = 0;
	}

	constraints_lowerbound[x_start]   = vars[x_start];
	constraints_lowerbound[y_start]   = vars[y_start];
	constraints_lowerbound[psi_start] = vars[psi_start];
	constraints_lowerbound[v_start]   = vars[v_start];
	constraints_lowerbound[cte_start] = vars[cte_start];
	constraints_lowerbound[epsi_start]= vars[epsi_start];

	constraints_upperbound[x_start]   = vars[x_start];
	constraints_upperbound[y_start]   = vars[y_start];
	constraints_upperbound[psi_start] = vars[psi_start];
	constraints_upperbound[v_start]   = vars[v_start];
	constraints_upperbound[cte_start] = vars[cte_start];
	constraints_upperbound[epsi_start]= vars[epsi_start];
	
	// object that computes objective and constraints
	FG_eval fg_eval(coeffs);

	//
	// NOTE: You don't have to worry about these options
	//
	// options for IPOPT solver
	std::string options;
	// Uncomment this if you'd like more print information
	options += "Integer print_level  0\n";
	// NOTE: Setting sparse to true allows the solver to take advantage
	// of sparse routines, this makes the computation MUCH FASTER. If you
	// can uncomment 1 of these and see if it makes a difference or not but
	// if you uncomment both the computation time should go up in orders of
	// magnitude.
	options += "Sparse  true        forward\n";
	options += "Sparse  true        reverse\n";
	// NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
	// Change this as you see fit.
	options += "Numeric max_cpu_time          0.5\n";

	// place to return solution
	CppAD::ipopt::solve_result<Dvector> solution;

	// solve the problem
	CppAD::ipopt::solve<Dvector, FG_eval>(
		options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
		constraints_upperbound, fg_eval, solution);

	// Check some of the solution values
	ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

	// Cost
	auto cost = solution.obj_value;
	std::cout << "Cost " << cost << std::endl;

	// TODO: Return the first actuator values. The variables can be accessed with
	// `solution.x[i]`.
	//
	// {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
	// creates a 2 element double vector.

	//std::cout << "here 5\n";

	for (int i = 1; i<N; i++) {
		glob_x[i - 1] = solution.x[i + x_start];
		glob_y[i - 1] = solution.x[i + y_start];
	}

	//std::cout << "here 6\n";
	for(int i; i<solution.x.size(); i++) { cout << "solution.x[" << i << "]   - " << solution.x[i] << endl; }
	//
	// Return the first actuator values. The variables can be accessed with
	return { solution.x[x_start     + 1], 
			 solution.x[y_start     + 1],
		     solution.x[psi_start   + 1], 
			 solution.x[v_start     + 1],
		     solution.x[cte_start   + 1], 
			 solution.x[epsi_start  + 1],
		     solution.x[delta_start    ], 
			 solution.x[a_start        ]};
}