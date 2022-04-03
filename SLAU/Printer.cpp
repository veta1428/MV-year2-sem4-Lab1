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
