#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <time.h>
#include <vector>

#include "alloc.hpp"

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double sigmaTum,
                  const double vascDens, const double sigmaVasc,
                  std::vector<bool> &inTum, std::vector<bool> &inVes);
void evalR(const int nMethod, const int nModel);
void model(const double *x, double *y, const int nrow,
           const int ncol, const int nlayer, const double cellSize,
           const std::vector<bool> &inTum, const std::vector<bool> &inVes);
void model(const double *x, double *y,  const int nrow,
           const int ncol, const int nlayer, const double cellSize,
           const std::vector<bool> &inTum, const std::vector<bool> &inVes,
           const std::string nFTumDens, const std::string nFTumVol,
           const std::string nFVascDens, const std::string nFKilledCells,
           const std::string nFCycle, const std::string nFHypDens,
           const std::string nFPO2Stat, const std::string nFVegfStat);
void morris(const int K, const int L, const int N, const int nOut,
            const int p, const double *x0, const double *h,
            double **mu, double **sigma, const std::string nFInTissueDim = "",
            const std::string nFInTum = "", const std::string nFInVes = "");
void morrisRT(const int N, const int p, const std::string nFRefParInt,
              const std::string nFInTissueDim, const std::string nFInTum,
              const std::string nFInVes);
void morrisToy(const int N, const int p, const std::string nFRefParInt);
void morrisVarRange(const int K, const int kp, const int L,
                    const int N, const int nOut, const int p,
                    const std::string nFRefParInt, double ***mu, double ***sigma,
                    const std::string nFInTissueDim = "", const std::string nFInTum = "",
                    const std::string nFInVes = "");
void morrisVarRangeRT(const int kp, const int L, const int N, const int p,
                      const std::string nFRefParInt, const std::string nFInTissueDim,
                      const std::string nFInTum, const std::string nFInVes);
void morrisVarRangeToy(const int kp, const int L, const int N, const int p,
                       const std::string nFRefParInt);
void readInFiles(const std::string nFInTissueDim, const std::string nFInTum,
                 const std::string nFInVes, int &nrow, int &ncol, int &nlayer,
                 double &cellSize, std::vector<bool> &inTum, std::vector<bool> &inVes);
void sobol(const int K, const int N, const int nOut,
           const double *x0, const double *h,
           double **SI, double **TSI,
           double ***SIConv, double ***TSIConv, const std::string nFInTissueDim = "",
           const std::string nFInTum = "", const std::string nFInVes = "");
void sobolFromFiles(int K);
void sobolRT(const int N, const std::string nFRefParInt, const std::string nFInTissueDim,
             const std::string nFInTum, const std::string nFInVes);
void sobolToy(const int N, const std::string nFRefParInt);
void toyModel(double *x, double *y);
void var1ParRange(const int kp, const int L, const std::string nRefParInt,
                  const std::string nFInTissueDim, const std::string nFInTum,
                  const std::string nFInVes);
void varErr(const std::string nFVarPar, const std::string nFMostRelPar,
            const std::string nFLeastPar, const std::string nFInTissueDim,
            const std::string nFInTum, const std::string nFInVes, const int L,
            const int P);
void varParFromFiles(const std::vector<std::string> nFPar, const std::string nFInTissueDim,
                     const std::string nFInTum, const std::string nFInVes);
void varStoch(const int N, const int P, const std::string nFRefParInt,
             const std::string nFInTissueDim, const std::string nFInTum,
             const std::string nFInVes);
