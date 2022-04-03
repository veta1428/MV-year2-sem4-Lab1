#ifndef _MATRIX_MANIP_REY_GUARD
#define _MATRIX_MANIP_REY_GUARD

#include "Models.h"

//TODO: may be can be optimized?
int ChangeLines(double** matrix, int first, int second, int rowLength);

int ChangeColumns(double** matrix, int first, int second, int columnLength);

double** MultMatrix(double** A, double** B, int n, int m, int k);

double* MultMatrixWithVector(double** matrix, double* vector, int n);

void DivideLineByElement(double** matrix, double element, int lineNumber, int lineLength, int startIndex);

//line element = line element - a * (elemnt of other line)
void LineSubstractAMOtherLine(double** matrix, int lineToSubstract, int lineSubstractFrom, double a, int rowLength, int startIndex);

int FindNotNullElementForGauss(double** matrix, int currentLine, int columnLength);

int FindBiggestInColumn(double** matrix, int lineStartFrom, int column, int columnLength);

int FindBiggestInRow(double** matrix, int columnStartFrom, int row, int rowLength);

double** CopyMatrix(double** A, int rows, int columns);

double CubeNorm(double** matrix, int lines);

double AbsSumLine(double** matrix, int line, int lineLength);

LUPFull SeparateLUP(LUP lup);

double Obuslovlennost(double** matrix, double** reversedMatrix, int n);

double Obuslovlennost(double** matrix, int n);

double* CopyVector(double* vector, int n);

#endif