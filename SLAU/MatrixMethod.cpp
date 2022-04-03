#include "MatrixMethod.h"
#include "MatrixManip.h"
#include <string.h>
#include "Models.h"
#include <cmath>
#include "Constants.h"

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
			LineSubstractAMOtherLine(matrix, i, j, mainElement, n, 0);

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