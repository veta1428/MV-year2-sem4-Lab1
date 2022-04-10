#ifndef _BENCHMARKS_REY_GUARD
#define _BENCHMARKS_REY_GUARD

#include "Models.h"
#include <string>
BenchmarkData ReportOne(double** matrix, int size, int variant, std::string filename);

BenchmarkData TestAll(int variant, int tests, int size);

ChangeBStatistics TestChangeB(double** matrix, int size, int variant);

void CheckWParam(std::string filename, double** matrix, int size, int variant);

#endif