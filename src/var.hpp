#ifndef DEF_SENSAN
#define DEF_SENSAN

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "alloc.hpp"
#include "evalModel.hpp"

void var1ParRange(const int kp, const int L, const std::string nRefParInt,
                  const std::string nFInTissueDim, const std::string nFInTum,
                  const std::string nFInVes, const std::string nFInPO2);
void varErr(const std::string nFVarPar, const std::string nFMostRelPar,
            const std::string nFLeastPar, const std::string nFInTissueDim,
            const std::string nFInTum, const std::string nFInVes,
            const std::string nFInPO2, const int L, const int P);
void varParFromFiles(const std::vector<std::string> nFPar,
                     const std::string nFInTissueDim,
                     const std::string nFInTum, const std::string nFInVes,
                     const std::string nFInPO2);
void varStoch(const int N, const int P, const std::string nFRefParInt,
              const std::string nFInTissueDim, const std::string nFInTum,
              const std::string nFInVes, const std::string nFInPO2);

#endif