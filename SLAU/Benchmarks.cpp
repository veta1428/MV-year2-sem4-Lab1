#include "Benchmarks.h"
#include <cmath>
#include "MatrixManip.h"
#include "Generator.h"
#include <string.h>
#include "MatrixMethod.h"
#include "Models.h"
#include <iostream>
#include <chrono>

void TestOne(double min, double max, int size)
{


}

BenchmarkData TestAll(int variant, int tests, int size) 
{
	double min = -pow(2, variant / 4);
	double max = -min;

	double obuslovlennostMIN = LONG_MAX;
	double obuslovlennostMAX = 0;
	double obuslovlennostAVERAGE = 0;
	double gaussDeltaMIN = LONG_MAX;
	double gaussDeltaMAX = 0;
	double gaussDeltaAVERAGE = 0;
	double lupDeltaMIN = LONG_MAX;
	double lupDeltaMAX = 0;
	double lupDeltaAVERAGE = 0;
	double ldltDeltaMIN = LONG_MAX;
	double ldltDeltaMAX = 0;
	double ldltDeltaAVERAGE = 0;
	double gaussTimeAVERAGE = 0;
	double buildLUPTimeAVERAGE = 0;
	double LUPTimeAVERAGE = 0;
	double LDLTTimeAVERAGE = 0;
	double relaxTimeAverage = 0;
	double iterMAX = 0;
	double iterMIN = LONG_MAX;
	double iterAVERAGE = 0;

	for (size_t i = 0; i < tests; i++)
	{
		double** testMatrix = AllocateMatrix(size, size);

		MakeRandomSimmetricMatrix(testMatrix, size, min, max);

		double* predefinedSolution = new double[size];
		memset(predefinedSolution, 0, sizeof(double) * size);

		MakeRandomVector(predefinedSolution, size, min, max);

		double* b = MultMatrixWithVector(testMatrix, predefinedSolution, size);


		double** testReversedMatrix = CopyMatrix(testMatrix, size, size);

		double** testGaussMatrix = CopyMatrix(testMatrix, size, size);
		double* testGaussVector = CopyVector(b, size);


		double** testLUPMatrix = CopyMatrix(testMatrix, size, size);
		double* testLUPVector = CopyVector(b, size);

		double** testLDLTMatrix = CopyMatrix(testMatrix, size, size);
		double* testLDLTVector = CopyVector(b, size);

		double** testRelaxMatrix = CopyMatrix(testMatrix, size, size);
		double* testRelaxVector = CopyVector(b, size);

		//reversed matrix
		double** reversedMatrix = GaussZordanReversedMatrix(testReversedMatrix, size);

		//obuslovlennost
		double obuslovlennost = Obuslovlennost(testMatrix, reversedMatrix, size);
		if (obuslovlennost < obuslovlennostMIN) {
			obuslovlennostMIN = obuslovlennost;
		}
		if (obuslovlennost > obuslovlennostMAX) {
			obuslovlennostMAX = obuslovlennost;
		}
		obuslovlennostAVERAGE += obuslovlennost;

		//gauss
		auto startGauss = std::chrono::steady_clock::now();
		double* gaussSolution = GaussLinearEq(testGaussMatrix, testGaussVector, size);
		auto endGauss = std::chrono::steady_clock::now();
		gaussTimeAVERAGE += std::chrono::duration_cast<std::chrono::nanoseconds>(endGauss - startGauss).count();

		double gaussDelta = CubeNormVector(Delta(gaussSolution, predefinedSolution, size), size);
		if (gaussDelta < gaussDeltaMIN) {
			gaussDeltaMIN = gaussDelta;
		}
		if (gaussDelta > gaussDeltaMAX) {
			gaussDeltaMAX = gaussDelta;
		}
		gaussDeltaAVERAGE += gaussDelta;

		//std::cout << "Gauss delta: " << gaussDelta << "\n";

		auto startBuildLUP= std::chrono::steady_clock::now();
		LUP lup = LUPByRow(testLUPMatrix, size);
		auto endBuildLUP = std::chrono::steady_clock::now();
		buildLUPTimeAVERAGE += std::chrono::duration_cast<std::chrono::nanoseconds>(endBuildLUP - startBuildLUP).count();


		auto startLUP = std::chrono::steady_clock::now();
		double* lupSolution = LUPSolveLinearEq(lup, testLUPVector);
		auto endLUP = std::chrono::steady_clock::now();
		LUPTimeAVERAGE += std::chrono::duration_cast<std::chrono::nanoseconds>(endLUP - startLUP).count();


		double lupDelta = CubeNormVector(Delta(lupSolution, predefinedSolution, size), size);
		if (lupDelta < lupDeltaMIN) {
			lupDeltaMIN = lupDelta;
		}
		if (lupDelta > lupDeltaMAX) {
			lupDeltaMAX = lupDelta;
		}
		lupDeltaAVERAGE += lupDelta;
		//std::cout << "LUP delta: " << lupDelta << "\n";

		LDL_T ldlt = LDLT(testLDLTMatrix, size);

		auto startLDLT = std::chrono::steady_clock::now();
		double* ldltSolution = LDLTLinearEq(ldlt, testLDLTVector);
		auto endLDLT = std::chrono::steady_clock::now();
		LDLTTimeAVERAGE += std::chrono::duration_cast<std::chrono::nanoseconds>(endLDLT - startLDLT).count();


		double ldltDelta = CubeNormVector(Delta(ldltSolution, predefinedSolution, size), size);
		if (ldltDelta < ldltDeltaMIN) {
			ldltDeltaMIN = ldltDelta;
		}
		if (ldltDelta > ldltDeltaMAX) {
			ldltDeltaMAX = ldltDelta;
		}
		ldltDeltaAVERAGE += ldltDelta;
		//std::cout << "LDLT delta: " << ldltDelta << "\n";

		auto startRelax = std::chrono::steady_clock::now();
		RelaxResult rr = RelaxIterations(testRelaxMatrix, testRelaxVector, size, 1 - (double)variant / 40);
		auto endRelax = std::chrono::steady_clock::now();
		relaxTimeAverage += std::chrono::duration_cast<std::chrono::nanoseconds>(endRelax - startRelax).count();

		double relaxDelta = CubeNormVector(rr.solution, size);
		int iter = rr.iterationAmount;

		if (iter < iterMIN) {
			iterMIN = iter;
		}
		if (iter > iterMAX) {
			iterMAX = iter;
		}
		iterAVERAGE += iter;
	}

	gaussTimeAVERAGE /= tests;
	buildLUPTimeAVERAGE /= tests;
	LUPTimeAVERAGE /= tests;
	LDLTTimeAVERAGE /= tests;
	obuslovlennostAVERAGE /= tests;
	gaussDeltaAVERAGE /= tests;
	lupDeltaAVERAGE /= tests;
	ldltDeltaAVERAGE /= tests;
	relaxTimeAverage /= tests;

	BenchmarkData bd;
	bd.buildLUPTimeAVERAGE = buildLUPTimeAVERAGE;
	bd.gaussTimeAVERAGE = gaussTimeAVERAGE;
	bd.LUPTimeAVERAGE = LUPTimeAVERAGE;
	bd.LDLTTimeAVERAGE = LDLTTimeAVERAGE;
	bd.obuslovlennostAVERAGE = obuslovlennostAVERAGE;
	bd.obuslovlennostMIN = obuslovlennostMIN;
	bd.obuslovlennostMAX = obuslovlennostMAX;
	bd.gaussDeltaAVERAGE = gaussDeltaAVERAGE;
	bd.lupDeltaAVERAGE = lupDeltaAVERAGE;
	bd.ldltDeltaAVERAGE = ldltDeltaAVERAGE;

	bd.gaussDeltaMAX = gaussDeltaMAX;
	bd.gaussDeltaMIN = gaussDeltaMIN;
	bd.ldltDeltaMAX = ldltDeltaMAX;
	bd.ldltDeltaMIN = ldltDeltaMIN;
	bd.lupDeltaMAX = lupDeltaMAX;
	bd.lupDeltaMIN = lupDeltaMIN;

	bd.iterAVERAGE = iterAVERAGE / tests;
	bd.iterMAX = iterMAX;
	bd.iterMIN = iterMIN;
	bd.relaxTimeAverage = relaxTimeAverage;

	return bd;

}