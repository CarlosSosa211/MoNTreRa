#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <time.h>
#include <vector>

#include "alloc.hpp"

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double sigmaTum,
                  const double vascDens, const double sigmaVasc,
                  std::vector<bool> &inTum, std::vector<bool> &inVes);
void evalR(const int nMethod, const int nModel);
void model(const double *x, double *y);
void model(const double *x, double *y, const std::string nFTumDens,
           const std::string nFTumVol, const std::string nFVascDens,
           const std::string nFKilledCells, const std::string nFCycle,
           const std::string nFHypDens, const std::string nFPO2Stat,
           const std::string nFVegfStat);
void morris(const int K, const int L, const int N, const int nOut,
            const int p, const double *x0, const double *h,
            double **mu, double **sigma);
void morrisRT(const int N, const int p, const std::string nFRefParInt);
void morrisToy(const int N, const int p, const std::string nFRefParInt);
void morrisVarRange(const int K, const int kp, const int L,
                    const int N, const int nOut, const int p,
                    const std::string nFRefParInt,
                    double ***mu, double ***sigma);
void morrisVarRangeRT(const int kp, const int L, const int N, const int p,
                      const std::string nFRefParInt);
void morrisVarRangeToy(const int kp, const int L, const int N, const int p,
                       const std::string nFRefParInt);
void sobol(const int K, const int N, const int nOut,
           const double *x0, const double *h,
           double **SI, double **TSI,
           double ***SIConv, double ***TSIConv);
void sobolFromFiles(int K);
void sobolRT(const int N, const std::string nFRefParInt);
void sobolToy(const int N, const std::string nFRefParInt);
void toyModel(double *x, double *y);
void var1Par(const int kp, const int p, const std::string nRefParMean);
