#include "Generator.h"
#include <random>

double fRand(double fMin, double fMax) 
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);

	//return rand() % 5;
}

void MakeRandomSimmetricMatrix(double** matrix, int n, double min, double max) {

	for (size_t i = 0; i < n; i++)
	{
		//because of the given formule
		matrix[i][i] = 1;
	}

	for (size_t i = 0; i < n; i++) //rows   i
	{
		for (size_t j = i + 1; j < n; j++) //columns   j
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
