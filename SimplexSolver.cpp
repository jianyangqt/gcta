/*
    Simple Simplex Solver Class
    Copyright (C) 2012  Tamas Bolner
	For more information, visit: http://blog.bolner.hu/2012/08/22/simplex/
	
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <Eigen/Dense>
#include "SimplexSolver.h"
#include <sstream>

using namespace Eigen;

/**
 * Constructor
 * 
 * @param int mode This can be one of these: SIMPLEX_MINIMIZE, SIMPLEX_MAXIMIZE
 * @param const VectorXf &objectiveFunction The coefficients of the objective function.
 * @param const MatrixXf &constraints Full matrix for the constraints. Contains also the righthand-side values.
 * @returns SimplexSolver
 */
SimplexSolver::SimplexSolver(int mode, const VectorXf &objectiveFunction, const MatrixXf &constraints) {
	int constantColumn, temp;
	this->foundSolution = false;
	this->optimum = 0;
	this->numberOfVariables = objectiveFunction.rows();
	
	/*
		Validate input parameters
	*/
	if (mode != SIMPLEX_MINIMIZE && mode != SIMPLEX_MAXIMIZE) {
		throw("SimplexSolver: invalid value for the mode parameter.");
	}

	if (objectiveFunction.rows() < 1) {
		throw("SimplexSolver: The coefficient vector of the objective function must contain at least one row.");
	}

	if (constraints.rows() < 1) {
		throw("SimplexSolver: The constraint matrix must contain at least one row.");
	}
	
	if (constraints.cols() != objectiveFunction.rows() + 1) {
        stringstream errmsg;
        errmsg<<"SimplexSolver: The constraint matrix has "<<constraints.cols()<<" columns, but should have "<<objectiveFunction.rows()+1<<", because the coefficient vector of the objective function has "<<objectiveFunction.rows()<<" rows."<<endl;
		throw(errmsg.str());
	}
	
	for (int i = 0; i < this->numberOfVariables; i++) {
		if (objectiveFunction(i) == 0) {
			throw("SimplexSolver: One of the coefficients of the objective function is zero.");
		}
	}

	temp = constraints.cols() - 1;
	for (int i = 0; i < constraints.rows(); i++) {
		if (constraints(i, temp) < 0) {
			throw("SimplexSolver: All righthand-side coefficients of the constraint matrix must be non-negative.");
		}
	}
	
	/*
		Build tableau
	*/
	if (mode == SIMPLEX_MAXIMIZE) {
		// Maximize
		this->tableau.resize(constraints.rows() + 1, this->numberOfVariables + constraints.rows() + 1);

		this->tableau <<	-objectiveFunction.transpose(),					MatrixXf::Zero(1, constraints.rows() + 1),
							constraints.leftCols(this->numberOfVariables),	MatrixXf::Identity(constraints.rows(), constraints.rows()), constraints.rightCols(1);
	} else {
		// Minimize: construct the Dual problem
		this->tableau.resize(this->numberOfVariables + 1, this->numberOfVariables + constraints.rows() + 1);

		this->tableau <<	-constraints.rightCols(1).transpose(),						MatrixXf::Zero(1, this->numberOfVariables + 1),
							constraints.leftCols(this->numberOfVariables).transpose(),	MatrixXf::Identity(this->numberOfVariables, this->numberOfVariables), objectiveFunction;
	}
	
	/*
		Simplex algorithm
	*/
	if (mode == SIMPLEX_MAXIMIZE) {
		// Maximize original problem
		if (!this->simplexAlgorithm(this->numberOfVariables)) {
			return;	// No solution
		}
	} else {
		// Maximize the dual of the minimization problem
		if (!this->simplexAlgorithm(constraints.rows())) {
			return;	// No solution
		}
	}
	
	/*
		Fetch solution
	*/
	constantColumn = this->tableau.cols() - 1;
	this->solution.resize(this->numberOfVariables);

	if (mode == SIMPLEX_MAXIMIZE) {
		// Maximize
		for (int i = 0; i < this->numberOfVariables; i++) {
			temp = this->getPivotRow(i);
			
			if (temp > 0) {
				// Basic variable
				this->solution(i) = this->tableau(temp, constantColumn);
			} else {
				// Non-basic variable
				this->solution(i) = 0;
			}
		}
		this->foundSolution = true;
		this->optimum = this->tableau(0, constantColumn);
	} else {
		// Minimize
		for (int i = 0; i < this->numberOfVariables; i++) {
			this->solution(i) = this->tableau(0, constraints.rows() + i);
		}
		this->foundSolution = true;
		this->optimum = this->tableau(0, constantColumn);
	}
}

/**
 * Returns true if a solution has been found.
 * Return false otherwise.
 *
 * @returns boolean
 */
bool SimplexSolver::hasSolution() {
	return this->foundSolution;
}

/**
 * Returns the maximum/minimum value of
 * the objective function.
 * 
 * @returns double
 */
double SimplexSolver::getOptimum() {
	return this->optimum;
}

/**
 * Returns a vector with the variable
 * values for the solution.
 *
 * return VectorXf
 */
VectorXf SimplexSolver::getSolution() {
	return this->solution;
}

/**
 * Searches for the pivot row in the given column, by calculating the ratios.
 * Tries to find smallest non-negative ratio.
 * Returns -1 if all possible pivots are 0 or if the ratios are negative.
 * Deals with cases like this:  0/negative < 0/positive
 * 
 * @param int column
 * @returns int Returns the number of the pivot row, or -1 if found none.
 */
int SimplexSolver::findPivot_min(int column) {
	int minIndex = -1;
	int constantColumn = this->tableau.cols() - 1;
	double minRatio = 0;
	double minConstant = 0;	// For the "0/negative < 0/positive" difference the constants have to be tracked also.
	double ratio;
	int rowNum = this->tableau.rows();
	
	for (int i = 1; i < rowNum; i++) {
		if (this->tableau(i, column) == 0) {
			continue;
		}

		ratio = this->tableau(i, constantColumn) / this->tableau(i, column);
		if (ratio < 0) {
			//The ratio must be non-negative
			continue;
		}

		if (minIndex == -1) {
			// First pivot candidate
			minIndex = i;
			minRatio = ratio;
			minConstant = this->tableau(i, constantColumn);
		} else {
			if (ratio == 0 && ratio == minRatio) {
				// 0/negative < 0/positive
				if (this->tableau(i, constantColumn) < minConstant) {
					minIndex = i;
					minRatio = ratio;
					minConstant = this->tableau(i, constantColumn);
				}
			} else if (ratio < minRatio) {
				minIndex = i;
				minRatio = ratio;
				minConstant = this->tableau(i, constantColumn);
			}
		}
	}

	return minIndex;
}

/**
 * Iterates through the this->tableau matrix to solve the problem.
 * 
 * @param int variableNum The number of variables (dimensions). (different for the minimization problem)
 * @returns bool Returns true if a solution has been found. Returns false otherwise.
 */
bool SimplexSolver::simplexAlgorithm(int variableNum) {
	MatrixXf::Index pivotColumn;
	int pivotRow;
	
	while (true) {
		/*
			Find pivot column, check for halt condition
		*/
		this->tableau.row(0).leftCols(variableNum).minCoeff(&pivotColumn);
		if (this->tableau(0, pivotColumn) >= 0) {
			//Found no negative coefficient
			break;
		}
		
		/*
			Find pivot row
		*/
		pivotRow = this->findPivot_min(pivotColumn);
		if (pivotRow == -1) {
			//no solution
			return false;
		}
		
		/*
			Do pivot operation
		*/
		this->tableau.row(pivotRow) /= this->tableau(pivotRow, pivotColumn);
		this->tableau(pivotRow, pivotColumn) = 1;	// For possible precision issues
		for (int i = 0; i < this->tableau.rows(); i++) {
			if (i == pivotRow) continue;
			
			this->tableau.row(i) -= this->tableau.row(pivotRow) * this->tableau(i, pivotColumn);
			this->tableau(i, pivotColumn) = 0;	// For possible precision issues
		}
	}

	return true;
}

/**
 * If the given column has only one coefficient with value 1 (except in topmost row), and all other
 * coefficients are zero, then returns the row of the non-zero value.
 * Otherwise return -1.
 * This method is used in the final step of maximization, when we read
 * the solution from the tableau.
 * 
 * @param int column
 * @returns int
 */
int SimplexSolver::getPivotRow(int column) {
	int one_row = -1;

	for (int i = 1; i < this->tableau.rows(); i++) {
		if (this->tableau(i, column) == 1) {
			if (one_row >= 0) {
				return -1;
			} else {
				one_row = i;
				continue;
			}
		} else if (this->tableau(i, column) != 0) {
			return -1;
		}
	}

	return one_row;
}
