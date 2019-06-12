#ifndef DEF_MORRIS
#define DEF_MORRIS

#include <random>
#include <string>
#include <vector>

#include "alloc.hpp"
#include "evalModel.hpp"

void morris(const int K, const int L, const int N, const int nOut, const int p,
            const double *x0, const double *h, double **mu, double **sigma,
            const std::string nFInTissueDim = "",
            const std::string nFInTum = "", const std::string nFInVes = "");
void morrisRT(const int N, const int p, const std::string nFRefParInt,
              const std::string nFInTissueDim, const std::string nFInTum,
              const std::string nFInVes);
void morrisToy(const int N, const int p, const std::string nFRefParInt);
void morrisVarRange(const int K, const int kp, const int L,
                    const int N, const int nOut, const int p,
                    const std::string nFRefParInt, double ***mu,
                    double ***sigma, const std::string nFInTissueDim = "",
                    const std::string nFInTum = "",
                    const std::string nFInVes = "");
void morrisVarRangeRT(const int kp, const int L, const int N, const int p,
                      const std::string nFRefParInt,
                      const std::string nFInTissueDim,
                      const std::string nFInTum, const std::string nFInVes);
void morrisVarRangeToy(const int kp, const int L, const int N, const int p,
                       const std::string nFRefParInt);

#endif
