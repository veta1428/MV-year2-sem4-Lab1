#ifndef _GENERATOR_REY_GUARD
#define _GENERATOR_REY_GUARD

//TODO: exclude fMax
double fRand(double fMin, double fMax);

void MakeRandomSimmetricMatrix(double** matrix, int n, double min, double max);

void MakeRandomVector(double* vector, int n, double min, double max);

#endif