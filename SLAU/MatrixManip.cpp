#include "MatrixManip.h"
#include <string.h>
#include "Constants.h"
#include <cmath>
#include "Models.h"
#include "MatrixMethod.h"

int ChangeLines(double** matrix, int first, int second, int rowLength)
{
	double temp = 0;
	for (size_t i = 0; i < rowLength; i++)
	{
		temp = matrix[first][i];
		matrix[first][i] = matrix[second][i];
		matrix[second][i] = temp;
	}
	return 0;
}

int ChangeColumns(double** matrix, int first, int second, int columnLength)
{
	double temp = 0;
	for (size_t i = 0; i < columnLength; i++)
	{
		temp = matrix[i][first];
		matrix[i][first] = matrix[i][second];
		matrix[i][second] = temp;
	}

	return 0;
}

double** MultMatrix(double** A, double** B, int n, int m, int k)
{
	double** C = new double* [n];

	for (size_t i = 0; i < n; i++)
	{
		C[i] = new double[n];
		memset(C[i], 0, n * sizeof(double));
	}

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < k; j++)
		{
			double sum = 0;

			for (size_t k = 0; k < m; k++)
			{
				sum += A[i][k] * B[k][j];
			}

			C[i][j] = sum;
		}
	}

	return C;
}

double* MultMatrixWithVector(double** matrix, double* vector, int n)
{
	double* b = new double[n];
	for (size_t i = 0; i < n; i++)
	{
		b[i] = 0;
		for (size_t j = 0; j < n; j++)
		{
			b[i] += matrix[i][j] * vector[j];
		}
	}
	return b;
}

void DivideLineByElement(double** matrix, double element, int lineNumber, int lineLength, int startIndex)
{
	for (size_t i = startIndex; i < lineLength; i++)
	{
		matrix[lineNumber][i] = matrix[lineNumber][i] / element;
	}
}

void LineSubstractAMOtherLine(double** matrix, int lineToSubstract, int lineSubstractFrom, double a, int n, int startIndex)
{
	for (size_t i = startIndex; i < n; i++)
	{
		matrix[lineSubstractFrom][i] -= a * matrix[lineToSubstract][i];
		double module = abs(matrix[lineSubstractFrom][i]);
		if (module < ZERO) {
			matrix[lineSubstractFrom][i] = 0;
		}
	}
}

int FindNotNullElementForGauss(double** matrix, int currentLine, int n)
{
	for (size_t i = currentLine + 1; i < n; i++)
	{
		if (matrix[i][currentLine] != 0)
		{
			return i;
		}
	}
	return -1;
}

int FindBiggestInColumn(double** matrix, int lineStartFrom, int column, int n)
{
	double biggest = 0;
	int lineNumber = 0;

	for (size_t i = lineStartFrom; i < n; i++)
	{
		if (abs(matrix[i][column]) > biggest)
		{
			biggest = abs(matrix[i][column]);
			lineNumber = i;
		}
	}

	return lineNumber;
}

int FindBiggestInRow(double** matrix, int columnStartFrom, int row, int n)
{
	double biggest = 0;
	int columnNumber = 0;

	for (size_t i = columnStartFrom; i < n; i++)
	{
		if (abs(matrix[row][i]) > biggest)
		{
			biggest = abs(matrix[row][i]);
			columnNumber = i;
		}
	}

	return columnNumber;
}

double** CopyMatrix(double** A, int rows, int columns)
{
	double** copy = new double* [rows];

	for (size_t i = 0; i < rows; i++)
	{
		copy[i] = new double[columns];
		for (size_t j = 0; j < columns; j++)
		{
			copy[i][j] = A[i][j];
		}
	}
	return copy;
}

LUPFull SeparateLUP(LUP lup)
{
	double** L = new double* [lup.n];
	double** U = new double* [lup.n];
	double** P = new double* [lup.n];

	for (size_t i = 0; i < lup.n; i++)
	{
		L[i] = new double[lup.n];
		U[i] = new double[lup.n];
		P[i] = new double[lup.n];

		memset(L[i], 0, sizeof(double) * lup.n);
		memset(U[i], 0, sizeof(double) * lup.n);
		memset(P[i], 0, sizeof(double) * lup.n);
	}

	//L
	for (size_t i = 0; i < lup.n; i++) //rows   i
	{
		for (size_t j = 0; j < lup.n; j++) //columns    j
		{
			double forPrint = 0;
			if (i > j)
			{
				forPrint = lup.LU[i][j];
			}
			else if (i == j) {
				forPrint = 1;
			}
			else {
				forPrint = 0;
			}
			L[i][j] = forPrint;
		}
	}

	//U
	for (size_t i = 0; i < lup.n; i++) //rows   i
	{
		for (size_t j = 0; j < lup.n; j++) //columns    j
		{
			double forPrint = 0;
			if (i <= j)
			{
				forPrint = lup.LU[i][j];
			}
			else {
				forPrint = 0;
			}
			U[i][j] = forPrint;
		}
	}

	//P
	for (size_t i = 0; i < lup.n; i++) //rows   i
	{
		int Erow = lup.P[i];
		for (size_t j = 0; j < lup.n; j++) //columns    j
		{
			double forPrint = 0;
			if (j == Erow)
			{
				forPrint = 1;
			}
			else {
				forPrint = 0;
			}
			P[i][j] = forPrint;
		}
	}

	LUPFull lf;
	lf.L = L;
	lf.U = U;
	lf.P = P;
	lf.n = lup.n;

	return lf;
}

double AbsSumLine(double** matrix, int line, int lineLength)
{
	double sum = 0;
	for (size_t i = 0; i < lineLength; i++)
	{
		sum += abs(matrix[line][i]);
	}
	return sum;
}

double CubeNorm(double** matrix, int lines)
{
	double maxSum = 0;
	for (size_t i = 0; i < lines; i++)
	{
		double lineSum = AbsSumLine(matrix, i, lines);

		if (lineSum > maxSum) {
			maxSum = lineSum;
		}
	}

	return maxSum;
}

double Obuslovlennost(double** matrix, double** reversedMatrix, int n)
{
	return CubeNorm(matrix, n) * CubeNorm(reversedMatrix, n);
}

double Obuslovlennost(double** matrix, int n)
{
	return Obuslovlennost(matrix, GaussZordanReversedMatrix(matrix, n), n);
}