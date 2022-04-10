#include "Printer.h"
#include <iostream>
#include <iomanip>

void PrintMatrix(double** matrix, int rows, int columns)
{
	for (size_t i = 0; i < rows; i++) //rows   i
	{
		std::cout << std::endl << "|";
		for (size_t j = 0; j < columns; j++) //columns    j
		{
			std::cout << std::setw(10) << matrix[i][j] << "    ";
		}
		std::cout << "|";
	}
}

void PrintMatrix(double** matrix, int rows)
{
	PrintMatrix(matrix, rows, rows);
}

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

void PrintBenchmarkResult(BenchmarkData bd) 
{
	std::cout << "\nWarning! If you see 1e+300 or -1e+300, that means No data for you here\n";
	std::cout << "\n\nTime (millisec):\n";
	std::cout << "Gauss: " << bd.gaussTimeAVERAGE << "\n";
	std::cout << "LUP build: " << bd.buildLUPTimeAVERAGE << "\n";
	std::cout << "LUP solve: " << bd.LUPTimeAVERAGE << "\n";
	std::cout << "LDLT solve: " << bd.LDLTTimeAVERAGE << "\n";
	std::cout << "Relax: " << bd.relaxTimeAverage << "\n";

	std::cout << "\n\nGauss:\n";
	std::cout << "Average delta: " << bd.gaussDeltaAVERAGE << "\n";
	std::cout << "Min delta: " << bd.gaussDeltaMIN << "\n";
	std::cout << "Max delta: " << bd.gaussDeltaMAX << "\n";

	std::cout << "\n\nLUP:\n";
	std::cout << "Average delta: " << bd.lupDeltaAVERAGE << "\n";
	std::cout << "Min delta: " << bd.lupDeltaMIN << "\n";
	std::cout << "Max delta: " << bd.lupDeltaMAX << "\n";

	std::cout << "\n\nLDLT\n";
	std::cout << "Average delta: " << bd.ldltDeltaAVERAGE << "\n";
	std::cout << "Min delta: " << bd.ldltDeltaMIN << "\n";
	std::cout << "Max delta: " << bd.ldltDeltaMAX << "\n";

	std::cout << "\n\nRelax\n";
	std::cout << "Average iter: " << bd.iterAVERAGE << "\n";
	std::cout << "Min iter: " << bd.iterMIN << "\n";
	std::cout << "Max iter: " << bd.iterMAX << "\n";
	std::cout << "Delta: " << bd.epsilonRelaxationsAVERAGE << "\n";

	std::cout << "\n\nObuslovlennost\n";
	std::cout << "Average ob: " << bd.obuslovlennostAVERAGE << "\n";
	std::cout << "Min ob: " << bd.obuslovlennostMIN << "\n";
	std::cout << "Max ob: " << bd.obuslovlennostMAX << "\n";
}
