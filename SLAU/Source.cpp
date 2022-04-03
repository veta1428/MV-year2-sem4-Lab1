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
	//float myVar = 11;
	//double min = -pow(2, myVar / 4);
	//double max = -min;

	//1
	double** matrix = new double*[n];
	for (size_t i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
		memset(matrix[i], 0, n * sizeof(*matrix[i]));
	}

	MakeRandomSimmetricMatrix(matrix, n, 0, 0);

	double** copyMatrixForReverse = CopyMatrix(matrix, n, n);

	double** matrixCopy = CopyMatrix(matrix, n, n);

	matrix[0][0] = 144;
	matrix[0][1] = -12;
	matrix[0][2] = -120;
	matrix[1][0] = -12;
	matrix[1][1] = -120;
	matrix[1][2] = -12;
	matrix[2][0] = -120;
	matrix[2][1] = -12;
	matrix[2][2] = 87;

	std::cout << "\nA =\n ";
	PrintMatrix(matrix, n);
	//double** matrixCopy = CopyMatrix(matrix, n, n);
	
	LDL_T ldlt = LDLT(matrix, n);

	std::cout << "\nprint ldlt\n";

	PrintMatrix(ldlt.LT, n);

	std::cout << "\n";
	for (size_t i = 0; i < n; i++)
	{
		std::cout << ldlt.isNegativeDiag[i] << " ";
	}
	
	double* solutionVector = new double[n];

	//solutionVector[0] = 1;
	//solutionVector[1] = -2;
	//solutionVector[2] = -1;
	MakeRandomVector(solutionVector, n, 0, 0);

	std::cout << "\nSolution vector\n";

	PrintVector(solutionVector, n);
	
	double* b = MultMatrixWithVector(matrix, solutionVector, n);
	double* bCopy = CopyVector(b, n);


	std::cout << "\n b vector\n";
	PrintVector(b, n);

	//LDL_T ldlt = LDLT(matrix, n);

	LUP lup = LUPByRow(matrix, n);

	double* solLUP = LUPSolveLinearEq(lup, b);
	std::cout << "\n Solution by LUP\n";

	PrintVector(solLUP, n);
	double* solGauss = GaussLinearEq(matrixCopy, bCopy, n);

	
	//PrintVector(solLUP, n);
	std::cout << "\n solution vector by gauss\n";

	PrintVector(solGauss, n);
	return 0;
}