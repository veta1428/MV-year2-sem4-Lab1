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

#endif