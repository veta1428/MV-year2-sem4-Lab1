#include "FileManager.h"
#include <string>
#include <fstream>
#include "MatrixManip.h"

void SaveMatrixToFile(std::string filename, int rows, int columns, double** matrix) 
{
	std::ofstream fout;
	fout.open(filename, std::ios::trunc);
	fout << rows << " " << columns << " ";
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < columns; j++)
		{
			fout << matrix[i][j] << " ";
		}
	}

	fout.close();
}

Matrix ReadMatrixFromFile(std::string file)
{
	std::ifstream fin;
	fin.open(file, std::ios::in);

	int row = 0;
	int column = 0;

	fin >> row >> column;

	double** matrix = AllocateMatrix(row, column);

	for (size_t i = 0; i < row; i++)
	{
		for (size_t j = 0; j < column; j++)
		{
			fin >> matrix[i][j];
		}
	}

	Matrix m;
	m.matrix = matrix;
	m.rows = row;
	m.columns = column;

	return m;
}