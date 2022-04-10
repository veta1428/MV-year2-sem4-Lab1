#include <iostream>
#include <iomanip>
#include "Benchmarks.h"
#include "Models.h"
#include "FileManager.h"
#include "Printer.h"
#include "MatrixManip.h"
#include "Generator.h"

#define MY_VARIANT 11
#define TEST_NUMBER 100
#define TEST_MATRIX_SIZE 255

int main() 
{
	srand(time(0));
	std::cout << "Testing. Please, do not interrupt...\n";
	BenchmarkData bd = TestAll(MY_VARIANT, TEST_NUMBER, TEST_MATRIX_SIZE);
	PrintBenchmarkResult(bd);

	Matrix a1 = ReadMatrixFromFile("a1.matrix");

	std::cout << "\n*******************A1 matrix test**************************\n";
	BenchmarkData a1b = ReportOne(a1.matrix, a1.rows, MY_VARIANT, "Tests\\A1_stats_");
	PrintBenchmarkResult(a1b);

	std::cout << "\n*******************A2 matrix test**************************\n";
	Matrix a2_ = ReadMatrixFromFile("a2.matrix");
	
	double** copy = CopyMatrix(a2_.matrix, a2_.rows, a2_.columns);
	double** a2_T = Transpone(copy, a2_.rows);
	
	double** a2 = MultMatrix(a2_T, a2_.matrix, a2_.rows, a2_.rows, a2_.rows);
	
	BenchmarkData a2b = ReportOne(a2, a2_.rows, MY_VARIANT, "Tests\\A2_stats_");
	PrintBenchmarkResult(a2b);
	
	std::cout << "\nTesting vector b change...\n";

	//std::cout << "\n*********************A2 b experiments************************\n";

	ChangeBStatistics s = TestChangeB(a2, a2_.rows, MY_VARIANT, "Tests\\A2");

	//PrintVector(s.solutionDelta, s.size);
	//PrintVector(s.bDelta, s.size);
	//std::cout << s.obuslovlennost;

	//std::cout << "\n*********************The worst matrix b experiments************************\n";

	Matrix worst = ReadMatrixFromFile("obuslovlennost.biggest");
	ChangeBStatistics sWorst = TestChangeB(worst.matrix, worst.rows, MY_VARIANT, "Tests\\Worst");

	//PrintVector(sWorst.solutionDelta, sWorst.size);
	//PrintVector(sWorst.bDelta, sWorst.size);
	//std::cout << sWorst.obuslovlennost;

	std::cout << "\nTesting relaxations method...\n";

	//relax 0.8 1.0 1.2
	CheckWParam("Tests\\A2", a2, a2_.rows, MY_VARIANT);
	Matrix worstRelax = ReadMatrixFromFile("obuslovlennost.biggest");

	CheckWParam("Tests\\Worst", worstRelax.matrix, worstRelax.rows, MY_VARIANT);
	std::cout << "\nDone. (results are in files)\n";
	return 0;
}