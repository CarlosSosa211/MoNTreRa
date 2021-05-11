#ifndef DEF_CREATEANDREADFILES
#define DEF_CREATEANDREADFILES

#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "simpcell.hpp"
#include "treatment.hpp"

static Treatment treat;

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double radRatioTum,
                  const double sigmaTum, const double vascDens,
                  const double sigmaVasc, std::vector<bool> &inTum,
                  std::vector<bool> &inVes);
int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double radRatioTum,
                  const double sigmaTum, const double vascDens,
                  std::vector<bool> &inTum, std::vector<bool> &inVes);
int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double vascDens, const std::vector<bool> &inTum,
                  std::vector<bool> &inVes);
int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double vascDens, const double sigmaVasc,
                  std::vector<bool> &inVes);
int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double vascDens,
                  std::vector<bool> &inTum, std::vector<bool> &inVes);
int createInFiles(const double cellSize, const double tumArea,
                  const double radRatioTum, const double tumDens,
                  const double vascDens, int &nrow, int &ncol, int &nlayer,
                  std::vector<bool> &inTum, std::vector<bool> &inVes);
void readInFiles(const std::string nFInTissueDim, const std::string nFInTum,
                 const std::string nFInVes, int &nrow, int &ncol, int &nlayer,
                 double &cellSize, std::vector<bool> &inTum,
                 std::vector<bool> &inVes, const std::string nFTreatment = "",
                 Treatment &Treatment = treat);
void readInFiles(const std::string nFInTissueDim, const std::string nFInTumVes,
                 int &nrow, int &ncol, int &nlayer, double &cellSize,
                 std::vector<bool> &inTumVes);
void readInFiles(const std::string nFInTissueDim, int &nrow, int &ncol,
                 int &nlayer, double &cellSize);
void readInFiles(const std::string nFInTissuePar, double &cellSize,
                 double &tumArea, double &radRatioTum, double &tumDens,
                 double &vascDens, const std::string nFTreatment = "",
                 Treatment &Treatment = treat);
void readInFilesOxy(const std::string nFInTissueOxy, bool &art, int &nrow,
                    int &ncol, int &nlayer, double &cellSize, double &vascDens,
                    double &sigmaVasc);
void readInFilesTCP(const std::string nFInTissueTCP,
                    const std::vector<std::string> nFTreatmentTCP, bool &art,
                    int &nrow, int &ncol, int &nlayer, double &cellSize,
                    double &tumDens, double &radRatioTum, double &sigmaTum,
                    double &vascDens, double &sigmaVasc,
                    std::vector<Treatment> &treatment);
void writeInFiles(const std::string nFInTissuePar);
void writeInVesFile(const std::string nFInTissueDim,
                    const std::string nFInTum);

#endif
