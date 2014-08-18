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

#pragma once
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#define SIMPLEX_MINIMIZE 1
#define SIMPLEX_MAXIMIZE 2

class SimplexSolver {
public:
	MatrixXf tableau;
	bool foundSolution;
	double optimum;
	VectorXf solution;
	int numberOfVariables;

	int findPivot_min(int column);
	bool simplexAlgorithm(int variableNum);
	int getPivotRow(int column);

	SimplexSolver(int mode, const VectorXf &objectiveFunction, const MatrixXf &constraints);
	bool hasSolution();
	double getOptimum();
	VectorXf getSolution();
};
