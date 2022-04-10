#include "Benchmarks.h"
#include <cmath>
#include "MatrixManip.h"
#include "Generator.h"
#include <string.h>
#include "MatrixMethod.h"
#include "Models.h"
#include <iostream>
#include <chrono>
#include "FileManager.h"
#include "Constants.h"
#include <fstream>

#define TEST_NUMBER 1

#define BIG_DOUBLE 1e+300
#define SMALL_DOUBLE -1e+300
#define CHANGE_B_DELTAS_AMOUNT 20

BenchmarkData ReportOne(double** testMatrix, int size, int variant, std::string filename)
{
	double min = -pow(2, variant / 4);
	double max = -min;

	double obuslovlennostMIN = BIG_DOUBLE;
	double obuslovlennostMAX = SMALL_DOUBLE;
	double obuslovlennostAVERAGE = 0;
	double gaussDeltaMIN = BIG_DOUBLE;
	double gaussDeltaMAX = SMALL_DOUBLE;
	double gaussDeltaAVERAGE = 0;
	double lupDeltaMIN = BIG_DOUBLE;
	double lupDeltaMAX = SMALL_DOUBLE;
	double lupDeltaAVERAGE = 0;
	double ldltDeltaMIN = BIG_DOUBLE;
	double ldltDeltaMAX = SMALL_DOUBLE;
	double ldltDeltaAVERAGE = 0;
	double gaussTimeAVERAGE = 0;
	double buildLUPTimeAVERAGE = 0;
	double LUPTimeAVERAGE = 0;
	double LDLTTimeAVERAGE = 0;
	double relaxTimeAverage = 0;
	double iterMAX = SMALL_DOUBLE;
	double iterMIN = BIG_DOUBLE;
	double iterAVERAGE = 0;
	double eAVERAGE = 0;
	for (size_t i = 0; i < TEST_NUMBER; i++)
	{

		double* predefinedSolution = new double[size];

		memset(predefinedSolution, 0, sizeof(double) * size);

		double* copyOfPredSol = CopyVector(predefinedSolution, size);

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
		gaussTimeAVERAGE += std::chrono::duration_cast<std::chrono::microseconds>(endGauss - startGauss).count();

		double gaussDelta = NormVector(Delta(gaussSolution, predefinedSolution, size), size);

		if (gaussDelta < gaussDeltaMIN) {
			gaussDeltaMIN = gaussDelta;
		}
		if (gaussDelta > gaussDeltaMAX) {
			gaussDeltaMAX = gaussDelta;
		}
		gaussDeltaAVERAGE += gaussDelta;

		auto startBuildLUP = std::chrono::steady_clock::now();
		LUP lup = LUPByRow(testLUPMatrix, size);
		auto endBuildLUP = std::chrono::steady_clock::now();
		buildLUPTimeAVERAGE += std::chrono::duration_cast<std::chrono::microseconds>(endBuildLUP - startBuildLUP).count();


		auto startLUP = std::chrono::steady_clock::now();
		double* lupSolution = LUPSolveLinearEq(lup, testLUPVector);
		auto endLUP = std::chrono::steady_clock::now();
		LUPTimeAVERAGE += std::chrono::duration_cast<std::chrono::microseconds>(endLUP - startLUP).count();


		double lupDelta = NormVector(Delta(lupSolution, predefinedSolution, size), size);
		if (lupDelta < lupDeltaMIN) {
			lupDeltaMIN = lupDelta;
		}
		if (lupDelta > lupDeltaMAX) {
			lupDeltaMAX = lupDelta;
		}
		lupDeltaAVERAGE += lupDelta;

		LDL_T ldlt = LDLT(testLDLTMatrix, size);

		auto startLDLT = std::chrono::steady_clock::now();
		double* ldltSolution = LDLTLinearEq(ldlt, testLDLTVector);
		auto endLDLT = std::chrono::steady_clock::now();
		LDLTTimeAVERAGE += std::chrono::duration_cast<std::chrono::microseconds>(endLDLT - startLDLT).count();


		double ldltDelta = NormVector(Delta(ldltSolution, predefinedSolution, size), size);
		if (ldltDelta < ldltDeltaMIN) {
			ldltDeltaMIN = ldltDelta;
		}
		if (ldltDelta > ldltDeltaMAX) {
			ldltDeltaMAX = ldltDelta;
		}
		ldltDeltaAVERAGE += ldltDelta;

		auto startRelax = std::chrono::steady_clock::now();
		RelaxResult rr = RelaxIterations(testRelaxMatrix, testRelaxVector, size, 1 - (double)variant / 40, true, copyOfPredSol, filename);
		auto endRelax = std::chrono::steady_clock::now();
		relaxTimeAverage += std::chrono::duration_cast<std::chrono::microseconds>(endRelax - startRelax).count();

		eAVERAGE += rr.e;
		int iter = rr.iterationAmount;

		if (iter < iterMIN) {
			iterMIN = iter;
		}
		if (iter > iterMAX) {
			iterMAX = iter;
		}
		iterAVERAGE += iter;
	}

	gaussTimeAVERAGE /= TEST_NUMBER;
	buildLUPTimeAVERAGE /= TEST_NUMBER;
	LUPTimeAVERAGE /= TEST_NUMBER;
	LDLTTimeAVERAGE /= TEST_NUMBER;
	obuslovlennostAVERAGE /= TEST_NUMBER;
	gaussDeltaAVERAGE /= TEST_NUMBER;
	lupDeltaAVERAGE /= TEST_NUMBER;
	ldltDeltaAVERAGE /= TEST_NUMBER;
	relaxTimeAverage /= TEST_NUMBER;
	eAVERAGE /= TEST_NUMBER;

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

	bd.iterAVERAGE = iterAVERAGE / TEST_NUMBER;
	bd.iterMAX = iterMAX;
	bd.iterMIN = iterMIN;
	bd.relaxTimeAverage = relaxTimeAverage;
	bd.epsilonRelaxationsAVERAGE = eAVERAGE;
	return bd;
}

BenchmarkData TestAll(int variant, int tests, int size) 
{
	double min = -pow(2, variant / 4);
	double max = -min;

	double obuslovlennostMIN = BIG_DOUBLE;
	double obuslovlennostMAX = SMALL_DOUBLE;
	double obuslovlennostAVERAGE = 0;
	double gaussDeltaMIN = BIG_DOUBLE;
	double gaussDeltaMAX = SMALL_DOUBLE;
	double gaussDeltaAVERAGE = 0;
	double lupDeltaMIN = BIG_DOUBLE;
	double lupDeltaMAX = SMALL_DOUBLE;
	double lupDeltaAVERAGE = 0;
	double ldltDeltaMIN = BIG_DOUBLE;
	double ldltDeltaMAX = SMALL_DOUBLE;
	double ldltDeltaAVERAGE = 0;
	double gaussTimeAVERAGE = 0;
	double buildLUPTimeAVERAGE = 0;
	double LUPTimeAVERAGE = 0;
	double LDLTTimeAVERAGE = 0;
	double relaxTimeAverage = 0;
	double iterMAX = SMALL_DOUBLE;
	double iterMIN = BIG_DOUBLE;
	double iterAVERAGE = 0;
	double eAVERAGE = 0;

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
			SaveMatrixToFile("obuslovlennost.biggest", size, size, testMatrix);
		}
		obuslovlennostAVERAGE += obuslovlennost;

		//gauss
		auto startGauss = std::chrono::steady_clock::now();
		double* gaussSolution = GaussLinearEq(testGaussMatrix, testGaussVector, size);
		auto endGauss = std::chrono::steady_clock::now();
		gaussTimeAVERAGE += std::chrono::duration_cast<std::chrono::microseconds>(endGauss - startGauss).count();

		double gaussDelta = NormVector(Delta(gaussSolution, predefinedSolution, size), size);
		if (gaussDelta < gaussDeltaMIN) {
			gaussDeltaMIN = gaussDelta;
		}
		if (gaussDelta > gaussDeltaMAX) {
			gaussDeltaMAX = gaussDelta;
		}
		gaussDeltaAVERAGE += gaussDelta;

		auto startBuildLUP= std::chrono::steady_clock::now();
		LUP lup = LUPByRow(testLUPMatrix, size);
		auto endBuildLUP = std::chrono::steady_clock::now();
		buildLUPTimeAVERAGE += std::chrono::duration_cast<std::chrono::microseconds>(endBuildLUP - startBuildLUP).count();


		auto startLUP = std::chrono::steady_clock::now();
		double* lupSolution = LUPSolveLinearEq(lup, testLUPVector);
		auto endLUP = std::chrono::steady_clock::now();
		LUPTimeAVERAGE += std::chrono::duration_cast<std::chrono::microseconds>(endLUP - startLUP).count();


		double lupDelta = NormVector(Delta(lupSolution, predefinedSolution, size), size);
		if (lupDelta < lupDeltaMIN) {
			lupDeltaMIN = lupDelta;
		}
		if (lupDelta > lupDeltaMAX) {
			lupDeltaMAX = lupDelta;
		}
		lupDeltaAVERAGE += lupDelta;

		LDL_T ldlt = LDLT(testLDLTMatrix, size);

		auto startLDLT = std::chrono::steady_clock::now();
		double* ldltSolution = LDLTLinearEq(ldlt, testLDLTVector);
		auto endLDLT = std::chrono::steady_clock::now();
		LDLTTimeAVERAGE += std::chrono::duration_cast<std::chrono::microseconds>(endLDLT - startLDLT).count();


		double ldltDelta = NormVector(Delta(ldltSolution, predefinedSolution, size), size);
		if (ldltDelta < ldltDeltaMIN) {
			ldltDeltaMIN = ldltDelta;
		}
		if (ldltDelta > ldltDeltaMAX) {
			ldltDeltaMAX = ldltDelta;
		}
		ldltDeltaAVERAGE += ldltDelta;

		double** L = AllocateMatrix(size, size);
		double** R = AllocateMatrix(size, size);
		auto startRelax = std::chrono::steady_clock::now();
		RelaxResult rr = RelaxIterations(testRelaxMatrix, testRelaxVector, size, 1 - (double)variant / 40);
		auto endRelax = std::chrono::steady_clock::now();
		relaxTimeAverage += std::chrono::duration_cast<std::chrono::microseconds>(endRelax - startRelax).count();
		eAVERAGE += rr.e;
		double relaxDelta = NormVector(rr.solution, size);
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
	eAVERAGE /= tests;

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
	bd.epsilonRelaxationsAVERAGE = eAVERAGE;
	return bd;

}

ChangeBStatistics TestChangeB(double** matrix, int size, int variant, std::string filename)
{
	std::ofstream foutDB;
	std::ofstream foutDS;

	std::string fncopy = std::string(filename);
	foutDB.open(filename.append(".deltaB"), std::ios::trunc);
	foutDS.open(fncopy.append(".deltaSolution"), std::ios::trunc);

	double* exactSolution = new double[size];
	double min = -pow(2, variant / 4);
	double max = -min;
	MakeRandomVector(exactSolution, size, min, max);

	double* b = MultMatrixWithVector(matrix, exactSolution, size);


	double* solutionDelta = new double[CHANGE_B_DELTAS_AMOUNT];
	double* bDelta = new double[CHANGE_B_DELTAS_AMOUNT];
	double** forObuslobl = CopyMatrix(matrix, size, size);

	double bNorm = NormVector(b, size);
	double solutionNorm = NormVector(exactSolution, size);

	double* solutionToCompare = CopyVector(exactSolution, size);
	double* bToCompare = CopyVector(b, size);

	double plus = 1.0 / CHANGE_B_DELTAS_AMOUNT;

	double exactSolutionNorm = NormVector(exactSolution, size);

	for (size_t i = 0; i < CHANGE_B_DELTAS_AMOUNT; i++)
	{
		ElementPlus(exactSolution, plus, size);
		double* newB = MultMatrixWithVector(matrix, exactSolution, size);

		double* deltaVectors = new double[size];
		MinusVectors(exactSolution, solutionToCompare, size, deltaVectors);

		double relDelta = NormVector(deltaVectors, size) / solutionNorm; /// bNorm;
		solutionDelta[i] = relDelta;

		foutDS << relDelta << "\n";

		double* deltaB = new double[size];
		MinusVectors(newB, bToCompare, size, deltaB);

		double bDeltaNorm = NormVector(deltaB, size) / bNorm;

		bDelta[i] = bDeltaNorm;
		foutDB << bDeltaNorm << "\n";
	}

	foutDB.close();
	foutDS.close();

	ChangeBStatistics s;
	s.bDelta = bDelta;
	s.solutionDelta = solutionDelta;
	s.size = CHANGE_B_DELTAS_AMOUNT;
	s.obuslovlennost = Obuslovlennost(forObuslobl, size);
	return s;
}

void CheckWParam(std::string filename, double** matrix, int size, int variant) 
{
	std::string f10 = std::string(filename);
	std::string f12 = std::string(filename);

	double* exactSolution = new double[size];
	double min = -pow(2, variant / 4);
	double max = -min;
	MakeRandomVector(exactSolution, size, min, max);

	double* b = MultMatrixWithVector(matrix, exactSolution, size);

	double** copy08 = CopyMatrix(matrix, size, size);
	double* b08 = CopyVector(b, size);

	double** copy10 = CopyMatrix(matrix, size, size);
	double* b10 = CopyVector(b, size);

	double** copy12 = CopyMatrix(matrix, size, size);
	double* b12 = CopyVector(b, size);

	RelaxIterations(copy08, b08, size, 0.8, true, exactSolution, filename.append("_0.8_"));
	RelaxIterations(copy10, b10, size, 1.0, true, exactSolution, f10.append("_1.0_"));
	RelaxIterations(copy12, b12, size, 1.2, true, exactSolution, f12.append("_1.2_"));
}