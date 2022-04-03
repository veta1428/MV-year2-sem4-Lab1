#ifndef _PRINTER_REY_GUARD
#define  _PRINTER_REY_GUARD

#include "Models.h"

void PrintMatrix(double** matrix, int rows, int columns);

void PrintMatrix(double** matrix, int rows);

//for debug
void PrintEquation(double** matrix, double* b, int n);

void PrintVector(double* b, int n);

void PrintLUP(LUP lup);

#endif