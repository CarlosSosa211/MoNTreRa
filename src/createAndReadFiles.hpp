#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "simpcell.hpp"
#include "treatment.hpp"

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double sigmaTum,
                  const double vascDens, const double sigmaVasc,
                  std::vector<bool> &inTum, std::vector<bool> &inVes);

void readInFiles(const std::string nFInTissueDim, const std::string nFInTum,
                 const std::string nFInVes, int &nrow, int &ncol, int &nlayer,
                 double &cellSize, std::vector<bool> &inTum, std::vector<bool> &inVes);

void readInFilesTCP(const std::string nFInTissueTCP,
                    const std::vector<std::string> nFTreatmentTCP,
                    int &nrow, int &ncol, int &nlayer, double &cellSize, double &tumDens,
                    double &sigmaTum, double &vascDens, double &sigmaVasc,
                    std::vector<Treatment> treatment);
