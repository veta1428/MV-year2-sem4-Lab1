#ifndef _MATRIX_METHOD_REY_GUARD
#define _MATRIX_METHOD_REY_GUARD

#include "Models.h"

double** GaussZordanReversedMatrix(double** matrix, int n);

double* GaussLinearEq(double** matrix, double* b, int n);

LUP LUPByRow(double** A, int n);

double* LUPSolveLinearEq(LUP lup, double* b);

LDL_T LDLT(double** matrix, int n);

double* LDLTLinearEq(LDL_T ldlt, double* b);

RelaxResult RelaxIterations(double** matrix, double* b, int size, double w);

#endif