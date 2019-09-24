#ifndef DEF_SOBOL
#define DEF_SOBOL

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "alloc.hpp"
#include "evalModel.hpp"

void sobol(const int K, const int N, const int nOut, const double *x0,
           const double *h, double **SI, double **TSI, double ***SIConv,
           double ***TSIConv, const std::string nFInTissueDim = "",
           const std::string nFInTum = "", const std::string nFInVes = "");
void sobolFromFiles(const int K);
void sobolRT(const int N, const std::string nFRefParInt,
             const std::string nFInTissueDim, const std::string nFInTum,
             const std::string nFInVes);
void sobolToy(const int N, const std::string nFRefParInt);

#endif
