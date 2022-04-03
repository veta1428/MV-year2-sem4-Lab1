#include <iostream>
#include <iomanip>

#define ZERO 0.00000000001

struct LUP 
{
	double** LU;
	int* P;
	int n;
};

struct LUPFull 
{
	double** L;
	double** U;
	double** P;
	int n;
};

//TODO: exclude fMax
double fRand(double fMin, double fMax)
{
	//double f = (double)rand() / RAND_MAX;
	//return fMin + f * (fMax - fMin);

	return rand() % 9;
}

double** MultMatrix(double** A, double** B, int n, int m, int k) 
{
	double** C = new double*[n];

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

void MakeSimmetricMatrix(double** matrix, int n, double min, double max) {

	for (size_t i = 0; i < n; i++)
	{
		matrix[i][i] = 1;//because of given formule
	}

	for (size_t i = 0; i < n; i++) //rows   i
	{
		for (size_t j = i + 1; j < n; j++) //columns    i
		{
			matrix[i][j] = fRand(min, max);
			matrix[i][i] += abs(matrix[i][j]);
			matrix[j][i] = matrix[i][j];
			matrix[j][j] += abs(matrix[j][i]);
		}
	}
}

void MakeRandomVector(double* vector, int n, double min, double max) {
	for (size_t i = 0; i < n; i++)
	{
		vector[i] = fRand(min, max);
	}
}

void PrintMatrix(double** matrix, int n) 
{
	for (size_t i = 0; i < n; i++) //rows   i
	{
		std::cout << std::endl << "|";
		for (size_t j = 0; j < n; j++) //columns    j
		{
			std::cout << std::setw(10) << matrix[i][j] << "    ";
		}
		std::cout << "|";
	}
}

//for debug
void PrintEquation(double** matrix, double* b, int n)
{
	for (size_t i = 0; i < n; i++) //rows   i
	{
		std::cout << std::endl << "|";
		for (size_t j = 0; j < n; j++) //columns    j
		{
			std::cout << std::setw(10) << matrix[i][j] << " ";
		}
		if (i == n / 2) {
			std::cout << "| * |" << std::setw(5) << "x" << i + 1 << " |";
			std::cout << "   =   |" << std::setw(10) << b[i] << " |";
		}
		else {
			std::cout << "|   |" << std::setw(5) << "x" << i + 1 << " |";
			std::cout << "       |" << std::setw(10) << b[i] << " |";
		}	
	}
}

void PrintVector(double* b, int n)
{
	std::cout << std::endl;
	for (size_t i = 0; i < n; i++)
	{
		std::cout << "|" << std::setw(10) << b[i] << "|" << std::endl;
	}
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

//TODO: may be can be optimized?
int ChangeLines(double** matrix, int first, int second, int n) 
{
	double temp = 0;
	for (size_t i = 0; i < n; i++)
	{
		temp = matrix[first][i];
		matrix[first][i] = matrix[second][i];
		matrix[second][i] = temp;
	}

	//double* temp = matrix[first];
	//matrix[first] = matrix[second];
	//matrix[second] = &temp;
	return 0;
}

int ChangeColumns(double** matrix, int first, int second, int n)
{
	double temp = 0;
	for (size_t i = 0; i < n; i++)
	{
		temp = matrix[i][first];
		matrix[i][first] = matrix[i][second];
		matrix[i][second] = temp;
	}

	return 0;
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

	//PrintMatrix(reversedMatrix, n);


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

		/*std::cout << "\nPrint matrix A after devision\n";

		PrintMatrix(matrix, n);*/
		DivideLineByElement(reversedMatrix, divideBy, i, n, 0);
		/*std::cout << "\nPrint matrix A-1 after devision\n";
		PrintMatrix(reversedMatrix, n);*/

		for (size_t j = i + 1; j < n; j++)
		{
			double mainElement = matrix[j][i];
			LineSubstractAMOtherLine(matrix, i, j, mainElement, n, 0);
			LineSubstractAMOtherLine(reversedMatrix, i, j, mainElement, n, 0);
		}

		/*std::cout << "\nPrint matrix A after minus\n";

		PrintMatrix(matrix, n);

		std::cout << "\nPrint matrix A-1 after minus\n";
		PrintMatrix(reversedMatrix, n);*/
	}

	/*std::cout << "\nA\n";

	PrintMatrix(matrix, n);

	std::cout << "\nA-1\n";
	PrintMatrix(reversedMatrix, n);*/
	//from bottom to top
	for (int i = n - 1; i >= 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			double mainElement = matrix[j][i];
			LineSubstractAMOtherLine(matrix, i, j, mainElement, n, 0);
			LineSubstractAMOtherLine(reversedMatrix, i, j, mainElement, n, 0);
		}

		/*std::cout << "\nPrint matrix A after minus\n";

		PrintMatrix(matrix, n);

		std::cout << "\nPrint matrix A-1 after minus\n";
		PrintMatrix(reversedMatrix, n);*/
	}

	return reversedMatrix;
}

int FindBiggestInColumn(double** matrix, int lineStartFrom, int column, int n) 
{
	double biggest = 0;
	int lineNumber = 0;

	for (size_t i = lineStartFrom; i < n; i++)
	{
		if (abs(matrix[i][column]) > biggest )
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


double* GaussLinearEq(double** matrix, double* b, int n) 
{
	//allocate memory for solution vector
	double* solution = new double[n];
	memset(solution, 0, sizeof(double) * n);

	//from top to bottom
	for (size_t i = 0; i < n; i++)
	{
		//prepare main element for actions
		//if (matrix[i][i] == 0)
		//{
		//	std::cout << "\nBefore find not zero\n";
		//	PrintMatrix(matrix, n);
		//	int changeWithLineIndex = FindNotNullElementForGauss(matrix, i, n);
		//	if (changeWithLineIndex == -1) {
		//		//TODO: so what?
		//		std::cout << "\nDeterminant is zero\n";
		//		return nullptr;
		//	}

		//	ChangeLines(matrix, i, changeWithLineIndex, n);
		//	
		//	double temp = b[i];
		//	b[i] = b[changeWithLineIndex];
		//	b[changeWithLineIndex] = temp;
		//}

		int lineToSwitchWith = FindBiggestInColumn(matrix, i, i, n);
		if (lineToSwitchWith != i) {
			ChangeLines(matrix, i, lineToSwitchWith, n);
		}

		if (matrix[i][i] == 0) {
			std::cout << "\nDeterminant is zero\n";
			return nullptr;
		}

		double divideBy = matrix[i][i];
		//divide line
		DivideLineByElement(matrix, divideBy, i, n, 0);

		std::cout << "\nPrint matrix A after devision\n";

		PrintMatrix(matrix, n);
		b[i] /= divideBy;
		/*std::cout << "\nPrint matrix A-1 after devision\n";
		PrintMatrix(reversedMatrix, n);*/

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

		std::cout << "\nPrint matrix A after minus\n";

		PrintMatrix(matrix, n);

		//std::cout << "\nPrint matrix A-1 after minus\n";
		//PrintMatrix(reversedMatrix, n);
	}
	std::cout << "\nTriangle\n";
	PrintMatrix(matrix, n);
	std::cout << "\n Vector b \n";
	PrintVector(b, n);

	for (int i = n - 1; i >= 0; i--)
	{
		int solved = n - 1 - i;

		double sumSolved = 0;

		for (int j = n - 1; j > i; j--)
		{
			//double sol = solution[j];
			//double matr = matrix[i][j];
			sumSolved += solution[j] * matrix[i][j];
		}
		//double newSol = b[i] - sumSolved;
		solution[i] = b[i] - sumSolved;
	}
	PrintVector(solution, n);
	/*std::cout << "\nA\n";

	PrintMatrix(matrix, n);

	std::cout << "\nA-1\n";
	PrintMatrix(reversedMatrix, n);*/
	//from bottom to top

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
			double a = A[j][i]/A[i][i];
			LineSubstractAMOtherLine(A, i, j, a, n, i);
			A[j][i] = a;
		}
	}

	PrintMatrix(A, n);
	LUP lup;
	lup.LU = A;
	lup.P = P;
	lup.n = n;
	for (size_t i = 0; i < n; i++)
	{
		std::cout << P[i] << " ";
	}
	return lup;
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

void PrintLUP(LUP lup)
{
	std::cout << "\n***********************LUP****************************\n";
	//print L
	std::cout << "\n\nL = \n\n";
	for (size_t i = 0; i < lup.n; i++) //rows   i
	{
		std::cout << std::endl << "|";
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
			std::cout << std::setw(10) << forPrint << "    ";
		}
		std::cout << "|";
	}

	//print U
	std::cout << "\n\nU = \n\n";
	for (size_t i = 0; i < lup.n; i++) //rows   i
	{
		std::cout << std::endl << "|";
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
			std::cout << std::setw(10) << forPrint << "    ";
		}
		std::cout << "|";
	}

	std::cout << "\n\nP = \n\n";
	for (size_t i = 0; i < lup.n; i++) //rows   i
	{
		std::cout << std::endl << "|";
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
			std::cout << std::setw(10) << forPrint << "    ";
		}
	}
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
	//PrintMatrix(GaussZordanReversedMatrix(matrix, n), n);

	//2

	std::cout << std::endl << "Solution =";
	//PrintVector(solutionVector, n);

	std::cout << std::endl << "Equation =";
	//PrintEquation(matrix, b, n);

	return 0;
}