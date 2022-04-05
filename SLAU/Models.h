#ifndef _MODELS_REY_GUARD
#define _MODELS_REY_GUARD

struct LUP
{
	double** LU;
	int* P;
	int n;
};

struct LUPFull
{
	double** L;
	double** U;
	double** P;
	int n;
};

struct LDL_T
{
	double** LT;
	bool* isNegativeDiag;
	int n;
};

struct BenchmarkData 
{
	double obuslovlennostMIN;
	double obuslovlennostMAX;
	double obuslovlennostAVERAGE;
	double gaussDeltaMIN;
	double gaussDeltaMAX;
	double gaussDeltaAVERAGE;
	double lupDeltaMIN;
	double lupDeltaMAX;
	double lupDeltaAVERAGE;
	double ldltDeltaMIN;
	double ldltDeltaMAX;
	double ldltDeltaAVERAGE;
	double gaussTimeAVERAGE;
	double buildLUPTimeAVERAGE;
	double LUPTimeAVERAGE;
	double LDLTTimeAVERAGE;
	double relaxTimeAverage;
	double iterMAX;
	double iterMIN;
	double iterAVERAGE;
};

struct ZeidelMatrix 
{
	double** L;
	double** R;
	double* g;
	int size;
};

struct CastZeidelToIteration 
{
	double** B_m;
	double* g_m;
};

struct RelaxResult {
	double* solution;
	int iterationAmount;
};

#endif