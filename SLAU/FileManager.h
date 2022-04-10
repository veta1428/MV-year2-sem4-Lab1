#ifndef _FILE_MANAGER_REY_GUARD
#define _FILE_MANAGER_REY_GUARD
#include <string>
#include "Models.h"
void SaveMatrixToFile(std::string filename, int rows, int columns, double** matrix);

Matrix ReadMatrixFromFile(std::string file);

#endif
