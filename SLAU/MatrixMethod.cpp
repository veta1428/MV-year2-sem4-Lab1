#include "MatrixMethod.h"
#include "MatrixManip.h"
#include <string.h>
#include "Models.h"
#include <cmath>
#include "Constants.h"
#include <iostream>

double** GaussZordanReversedMatrix(double** matrix, int n)
{
	//allocate reversed matrix
	double** reversedMatrix = new double* [n];
	for (size_t i = 0; i < n; i++)
	{
		reversedMatrix[i] = new double[n];
		memset(reversedMatrix[i], 0, sizeof(double) * n);
		reversedMatrix[i][i] = 1;
	}

	//from top to bottom
	for (size_t i = 0; i < n; i++)
	{
		//prepare main element for actions
		if (matrix[i][i] == 0)
		{
			int changeWithLineIndex = FindNotNullElementForGauss(matrix, i, n);
			if (changeWithLineIndex == -1) {
				return nullptr;
			}
			ChangeLines(matrix, i, changeWithLineIndex, n);
			ChangeLines(reversedMatrix, i, changeWithLineIndex, n);
		}

		double divideBy = matrix[i][i];

		//divide line
		DivideLineByElement(matrix, divideBy, i, n, 0);
		DivideLineByElement(reversedMatrix, divideBy, i, n, 0);

		for (size_t j = i + 1; j < n; j++)
		{
			double mainElement = matrix[j][i];
			LineSubstractAMOtherLine(matrix, i, j, mainElement, n, 0);
			LineSubstractAMOtherLine(reversedMatrix, i, j, mainElement, n, 0);
		}
	}

	for (int i = n - 1; i >= 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			double mainElement = matrix[j][i];
			LineSubstractAMOtherLine(matrix, i, j, mainElement, n, 0);
			LineSubstractAMOtherLine(reversedMatrix, i, j, mainElement, n, 0);
		}
	}

	return reversedMatrix;
}

double* GaussLinearEq(double** matrix, double* b, int n)
{
		//allocate memory for solution vector
	double* solution = new double[n];
	memset(solution, 0, sizeof(double) * n);

	//from top to bottom
	for (size_t i = 0; i < n; i++)
	{
		int lineToSwitchWith = FindBiggestInColumn(matrix, i, i, n);
		if (lineToSwitchWith != i) {
			ChangeLines(matrix, i, lineToSwitchWith, n);
		}

		if (matrix[i][i] == 0) {
			return nullptr;
		}

		double divideBy = matrix[i][i];
		//divide line
		DivideLineByElement(matrix, divideBy, i, n, 0);

		b[i] /= divideBy;


		for (size_t j = i + 1; j < n; j++)
		{
			double mainElement = matrix[j][i];
			LineSubstractAMOtherLine(matrix, i, j, mainElement, n, i);

			b[j] -= mainElement * b[i];

			//TODO: do we need to?
			if (abs(b[j]) < ZERO) 
			{
				b[j] = 0;
			}
		}
	}


	for (int i = n - 1; i >= 0; i--)
	{
		int solved = n - 1 - i;

		double sumSolved = 0;

		for (int j = n - 1; j > i; j--)
		{
			sumSolved += solution[j] * matrix[i][j];
		}

		solution[i] = b[i] - sumSolved;

		if (abs(solution[i]) < ZERO)
		{
			solution[i] = 0;
		}
	}

	return solution;
}

LUP LUPByRow(double** A, int n)
{
	int* P = new int[n];
	for (size_t i = 0; i < n; i++)
	{
		P[i] = i;
	}


	for (size_t i = 0; i < n; i++)
	{
		//change columns if needed
		int changeWithColumn = FindBiggestInRow(A, i, i, n);
		if (i != changeWithColumn) {
			ChangeColumns(A, i, changeWithColumn, n);

			int temp = P[i];
			P[i] = P[changeWithColumn];
			P[changeWithColumn] = temp;
		}

		for (size_t j = i + 1; j < n; j++)
		{
			double a = A[j][i] / A[i][i];
			LineSubstractAMOtherLine(A, i, j, a, n, i);
			A[j][i] = a;
		}
	}

	LUP lup;
	lup.LU = A;
	lup.P = P;
	lup.n = n;

	return lup;
}

double* LUPSolveLinearEq(LUP lup, double* b) 
{
	//LU * x' = b

	//Ly = b

	double* ySolution = new double[lup.n];
	memset(ySolution, 0, sizeof(double) * lup.n);


	for (int i = 0; i < lup.n; i++)
	{
		int solved = i;

		double sumSolved = 0;

		for (int j = 0; j < solved; j++)
		{
			sumSolved += ySolution[j] * lup.LU[i][j];
		}

		ySolution[i] = b[i] - sumSolved;
	}

	// U x' = ySolution
	//

	double* xMSolution = new double[lup.n];
	memset(xMSolution, 0, sizeof(double) * lup.n);

	for (int i = lup.n - 1; i >= 0; i--)
	{
		int solved = lup.n - 1 - i;

		double sumSolved = 0;

		for (int j = lup.n - 1; j > i; j--)
		{
			sumSolved += xMSolution[j] * lup.LU[i][j];
		}

		xMSolution[i] = (ySolution[i] - sumSolved) / lup.LU[i][i];
	}

	return xMSolution;
}

LDL_T LDLT(double** matrix, int n) 
{
	bool* isNegative = new bool[n];
	memset(isNegative, 0, sizeof(bool) * n);

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = i + 1; j < n; j++)
		{
			double a = matrix[j][i] / matrix[i][i];
			LineSubstractAMOtherLine(matrix, i, j, a, n, i);
		}
	}

	//divide
	for (size_t i = 0; i < n; i++)
	{
		double diag = matrix[i][i];
		double sqrtv = 1;
		if (diag < 0) {
			isNegative[i] = true;
			diag *= -1;
			sqrtv *= -1;
		}

		sqrtv *= sqrt(diag);

		for (size_t j = i; j < n; j++)
		{
			matrix[i][j] /= sqrtv;
		}
	}

	LDL_T t;
	t.n = n;
	t.isNegativeDiag = isNegative;
	t.LT = matrix;

	return t;
}

double* LDLTLinearEq(LDL_T ldlt, double* b) 
{
	// LD * Ltx = b
	// Ltx = y;
	// y' = Dy
	// L * y' = b
	double* yMSolution = new double[ldlt.n];
	memset(yMSolution, 0, sizeof(double) * ldlt.n);

	for (int i = 0; i < ldlt.n; i++)
	{
		int solved = i;

		double sumSolved = 0;

		for (int j = 0; j < solved; j++)
		{
			double sol = yMSolution[j];
			double ldly = ldlt.LT[j][i];
			sumSolved += yMSolution[j] * ldlt.LT[j][i];
		}

		double s = (b[i] - sumSolved) / ldlt.LT[i][i];
		yMSolution[i] = (b[i] - sumSolved) / ldlt.LT[i][i];
	}


	// y' = Dy
	double* ySolution = new double[ldlt.n];
	memset(ySolution, 0, sizeof(double) * ldlt.n);

	for (size_t i = 0; i < ldlt.n; i++)
	{
		ySolution[i] = yMSolution[i];
		if (ldlt.isNegativeDiag[i] == true) {
			ySolution[i] *= -1;
		}
	}

	// Ltx = y;
	double* xSolution = new double[ldlt.n];
	memset(xSolution, 0, sizeof(double) * ldlt.n);

	for (int i = ldlt.n - 1; i >= 0; i--)
	{
		int solved = ldlt.n - 1 - i;

		double sumSolved = 0;

		for (int j = ldlt.n - 1; j > i; j--)
		{
			double sol = xSolution[j];
			double ldly = ldlt.LT[i][j];
			sumSolved += xSolution[j] * ldlt.LT[i][j];
		}

		double s = (ySolution[i] - sumSolved) / ldlt.LT[i][i];
		xSolution[i] = (ySolution[i] - sumSolved) / ldlt.LT[i][i];
	}

	return xSolution;
}

RelaxResult RelaxIterations(double** matrix, double* b, int size, double w) 
{
	ZeidelMatrix zm = GetZeidelMatrixFromStrongDiagMatrix(matrix, size, b);
	CastZeidelToIteration cast = GetIterationMatrixFromRelax(zm, w);

	double* solution = new double[size];
	memset(solution, 1, sizeof(double) * size);

	int counter = 0;
	while (true) 
	{
		counter++;
		double* multRes = MultMatrixWithVector(cast.B_m, solution, size);

		SumVectors(multRes, cast.g_m, size);

		MinusVectors(solution, multRes, size);

		double norm = CubeNormVector(solution, size);

		solution = multRes;
		if (norm < E) 
		{
			break;
		}
	}

	RelaxResult rr;

	rr.solution = solution;
	rr.iterationAmount = counter;

	return rr;
}