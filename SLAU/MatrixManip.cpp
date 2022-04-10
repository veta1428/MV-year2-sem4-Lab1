#include "MatrixManip.h"
#include <string.h>
#include "Constants.h"
#include <cmath>
#include "Models.h"
#include "MatrixMethod.h"
#include <iostream>
#include "Printer.h"

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

double Norm(double** matrix, int lines)
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
	//std::cout << "\n" << Norm(matrix, n) << "\n" << Norm(reversedMatrix, n) << "\n";
	return Norm(matrix, n) * Norm(reversedMatrix, n);
}

double Obuslovlennost(double** matrix, int n)
{
	double** copy = CopyMatrix(matrix, n, n);
	return Obuslovlennost(matrix, GaussZordanReversedMatrix(copy, n), n);
}

double* CopyVector(double* vector, int n) 
{
	double* copy = new double[n];
	memset(copy, 0, sizeof(double) * n);

	for (size_t i = 0; i < n; i++)
	{
		copy[i] = vector[i];
	}
	return copy;
}

double** AllocateMatrix(int rows, int columns) {
	double** matrix = new double* [rows];

	for (size_t i = 0; i < columns; i++)
	{
		matrix[i] = new double[columns];
		memset(matrix[i], 0, sizeof(double) * columns);
	}

	return matrix;
}

double* Delta(double* firstVector, double* secondVector, int size) 
{
	double* delta = new double[size];
	memset(delta, 0, sizeof(double) * size);

	for (size_t i = 0; i < size; i++)
	{
		delta[i] = firstVector[i] - secondVector[i];
	}
	return delta;
}

double NormVector(double* vector, int size) 
{
	double norm = 0;

	for (size_t i = 0; i < size; i++)
	{
		if (abs(vector[i]) > norm)
		{
			norm = abs(vector[i]);
		}
	}

	return norm;
}

ZeidelMatrix GetZeidelMatrixFromStrongDiagMatrix(double** matrix, int size, double* b, double** L, double** R)
{
	for (size_t i = 0; i < size; i++)
	{
		b[i] = b[i] / matrix[i][i];
		for (size_t j = 0; j < size; j++)
		{
			if (i < j)
			{
				R[i][j] = -matrix[i][j] / matrix[i][i];
			}
			else if (i > j) {
				L[i][j] = -matrix[i][j] / matrix[i][i];
			}
		}
	}

	ZeidelMatrix zm;
	zm.L = L;
	zm.R = R;
	zm.size = size;
	zm.g = b;

	return zm;
}

CastZeidelToIteration GetIterationMatrixFromRelax(ZeidelMatrix zm, double w) 
{
	//xk+1 = (E-wL)^-1 * ((wR + E - wE)xk + wg)
	//       |_______|   |___________|    |__|
	//           L_            R_          


	//L_
	for (size_t i = 0; i < zm.size; i++)
	{
		zm.L[i][i] = 1;
		zm.g[i] *= w;
		zm.R[i][i] = 1 - w;
		for (size_t j = 0; j < zm.size; j++)
		{
			if (i > j)
			{
				zm.L[i][j] = -zm.L[i][j] * w;
			} else if (i < j)
			{
				zm.R[i][j] = zm.R[i][j] * w;
			}
		}
	}

	double** L_ = GaussZordanReversedMatrix(zm.L, zm.size);

	////R_, g_
	//for (size_t i = 0; i < zm.size; i++)
	//{
	//	zm.g[i] *= w;
	//	zm.R[i][i] = 1 - w;
	//	for (size_t j = 0; j < zm.size; j++)
	//	{
	//		if (i < j)
	//		{
	//			zm.R[i][j] = zm.R[i][j] * w;
	//		}
	//	}
	//}

	double** B_m = MultMatrix(L_, zm.R, zm.size, zm.size, zm.size);

	double* g_m = MultMatrixWithVector(L_, zm.g, zm.size);
	
	CastZeidelToIteration cast;
	cast.B_m = B_m;
	cast.g_m = g_m;

	return cast;
}
void SumVectors(double* destination, double* source, int size) 
{
	for (size_t i = 0; i < size; i++)
	{
		destination[i] += source[i];
	}
}


void MinusVectors(double* first, double* second, int size, double* result)
{
	for (size_t i = 0; i < size; i++)
	{
		result[i] = second[i] - first[i];
	}
}

double** Transpone(double** matrix, int size) 
{
	double** transponed = AllocateMatrix(size, size);

	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			transponed[i][j] = matrix[j][i];
		}
	}

	return transponed;
}

void ElementPlus(double* vector, double plus, int size) 
{
	for (size_t i = 0; i < size; i++)
	{
		vector[i] += plus;
	}
}