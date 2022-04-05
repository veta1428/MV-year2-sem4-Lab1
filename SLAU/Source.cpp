#include <iostream>
#include <iomanip>
#include "Benchmarks.h"
#include "Models.h"

int main() 
{
	BenchmarkData bd = TestAll(11, 5, 255);
	std::cout << "Time:\n";
	std::cout << "Gauss: " << bd.gaussTimeAVERAGE << "\n";
	std::cout << "LUP build: " << bd.buildLUPTimeAVERAGE << "\n";
	std::cout << "LUP solve: " << bd.LUPTimeAVERAGE << "\n";
	std::cout << "LDLT solve: " << bd.LDLTTimeAVERAGE << "\n";
	std::cout << "Relax: " << bd.relaxTimeAverage << "\n";

	std::cout << "\n\nGauss:\n";
	std::cout << "Average delta: " << bd.gaussDeltaAVERAGE << "\n";
	std::cout << "Min delta: " << bd.gaussDeltaMIN << "\n";
	std::cout << "Max delta: " << bd.gaussDeltaMAX << "\n";

	std::cout << "\n\nLUP:\n";
	std::cout << "Average delta: " << bd.lupDeltaAVERAGE << "\n";
	std::cout << "Min delta: " << bd.lupDeltaMIN << "\n";
	std::cout << "Max delta: " << bd.lupDeltaMAX << "\n";

	std::cout << "\n\nLDLT\n";
	std::cout << "Average delta: " << bd.ldltDeltaAVERAGE << "\n";
	std::cout << "Min delta: " << bd.ldltDeltaMIN << "\n";
	std::cout << "Max delta: " << bd.ldltDeltaMAX << "\n";

	std::cout << "\n\nRelax\n";
	std::cout << "Average iter: " << bd.iterAVERAGE << "\n";
	std::cout << "Min iter: " << bd.iterMIN << "\n";
	std::cout << "Max iter: " << bd.iterMAX << "\n";

	std::cout << "\n\nObuslovlennost\n";
	std::cout << "Average ob: " << bd.obuslovlennostAVERAGE<< "\n";
	std::cout << "Min ob: " << bd.obuslovlennostMIN << "\n";
	std::cout << "Max ob: " << bd.obuslovlennostMAX << "\n";
	return 0;
}