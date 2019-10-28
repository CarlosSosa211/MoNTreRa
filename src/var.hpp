#ifndef DEF_SENSAN
#define DEF_SENSAN

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

#include "alloc.hpp"
#include "evalModel.hpp"

void var1ParRange(const int kp, const int L, const std::string nRefParInt,
                  const std::string nFInTissueDim, const std::string nFInTum,
                  const std::string nFInVes);
void varArtTissue(const int P, const std::string nFArt,
                  const std::string nFInTissueDim,
                  const std::string nFRefParMean);
void varArtTissue(const int N, const int P, const std::string nFDensInt,
                  const std::string nFInTissueDim,
                  const std::string nFRefParMean);
void varArtTissue(const int N, const int P, const std::string nFDensInt,
                  const std::string nFSigmaTumInt,
                  const std::string nFInTissueDim,
                  const std::string nFRefParMean);
void varErr(const std::string nFVarPar, const std::string nFMostRelPar,
            const std::string nFLeastPar, const std::string nFInTissueDim,
            const std::string nFInTum, const std::string nFInVes, const int L,
            const int P);
void varParFromFiles(const std::vector<std::string> nFPar,
                     const std::string nFInTissueDim,
                     const std::string nFInTum, const std::string nFInVes);
void varStoch(const int N, const int P, const std::string nFRefParInt,
              const std::string nFInTissueDim, const std::string nFInTum,
              const std::string nFInVes);

#endif
