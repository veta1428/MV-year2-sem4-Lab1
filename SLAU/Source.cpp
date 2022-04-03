#include <iostream>
#include <iomanip>
#include "Models.h"
#include "Constants.h"
#include "Generator.h"
#include "Printer.h"
#include "MatrixManip.h"
#include "MatrixMethod.h"

int main() 
{
	int n = 3;
	float myVar = 11;
	double min = -pow(2, myVar / 4);
	double max = -min;

	//1
	double** matrix = new double*[n];
	for (size_t i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
		memset(matrix[i], 0, n * sizeof(*matrix[i]));
	}

	//MakeSimmetricMatrix(matrix, n, min, max);
	matrix[0][0] = 1;
	matrix[0][1] = 2;
	matrix[0][2] = -10;
	matrix[1][0] = 3;
	matrix[1][1] = -1;
	matrix[1][2] = 50;
	matrix[2][0] = 7;
	matrix[2][1] = 3;
	matrix[2][2] = 1;

	std::cout << "A =";
	PrintMatrix(matrix, n);

	double* solutionVector = new double[n];
	MakeRandomVector(solutionVector, n, min, max);
	//std::cout << "\nSolution\n";
	//PrintVector(solutionVector, n);
	LUP lup = LUPByRow(matrix, n);
	LUPFull lupFull = SeparateLUP(lup);
	std::cout << "\nChecked L * U\n";

	PrintMatrix(MultMatrix(lupFull.L, lupFull.U, n, n, n), n);
	double* b = MultMatrixWithVector(matrix, solutionVector, n);
	double* sol = GaussLinearEq(matrix, b, n);
	PrintVector(sol, n);
	PrintVector(solutionVector, n);
	return 0;
}