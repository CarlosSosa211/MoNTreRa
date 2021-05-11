#ifndef DEF_EVALMODEL
#define DEF_EVALMODEL

#include "constOxyTissue.hpp"
#include "createAndReadFiles.hpp"
#include "oxyTissue.hpp"
#include "rootSimulator.hpp"
#include "simpcell.hpp"
#include "tissue.hpp"
#include "treatment.hpp"


void model(const double *x, double *y, const int nrow,
           const int ncol, const int nlayer, const double cellSize,
           const std::vector<bool> &inTum, const std::vector<bool> &inVes);
void model(const double *x, double *y, const int nrow,
           const int ncol, const int nlayer, const double cellSize);
void model(const double *x, const int nrow, const int ncol, const int nlayer,
           const double cellSize, const std::vector<bool> &inVes,
           const std::string nFPO2);
void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize,
           const std::vector<bool> &inTum, const std::vector<bool> &inVes,
           const std::string nFTumDens, const std::string nFTumVol,
           const std::string nFVascDens, const std::string nFKilledCells,
           const std::string nFDeadDens, const std::string nFCycle,
           const std::string nFHypDens, const std::string nFPO2Stat,
           const std::string nFVegfStat);
void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize,
           const std::vector<bool> &inTum, const std::vector<bool> &inVes,
           Treatment *const treatment, const std::string nFTumDens,
           const std::string nFTumVol, const std::string nFVascDens,
           const std::string nFKilledCells, const std::string nFDeadDens,
           const std::string nFCycle, const std::string nFHypDens,
           const std::string nFPO2Stat, const std::string nFVegfStat,
           const std::string nFState = "", const std::string nFTimer = "",
           const std::string nFPO2 = "", const std::string nFVegf = "");
void oxy(const int N, const std::string nFInTissueOxy,
         const std::string nFParOxy, const std::string nFInTissueDim = "",
         const std::string nFInVes = "");
void reducedFracModel(const double *x, double *y, const int nrow,
                  const int ncol, const int nlayer, const double cellSize,
                  const std::vector<bool> &inTum,
                  const std::vector<bool> &inVes);
void reducedModel(const double *x, double *y, const int nrow, const int ncol,
                  const int nlayer, const double cellSize,
                  const std::vector<bool> &inTum,
                  const std::vector<bool> &inVes);
void reducedModel(const double *x, double *y, const int nrow, const int ncol,
                  const int nlayer, const double cellSize,
                  const std::vector<bool> &inTum,
                  const std::vector<bool> &inVes, Treatment *const treatment);
void reducedArtModel(const double *x, double *y);
void reducedArtModel(const double *x, double *y, const std::string nFTumDens,
                     const std::string nFTumVol, const std::string nFVascDens,
                     const std::string nFKilledCells,
                     const std::string nFDeadDens, const std::string nFCycle,
                     const std::string nFHypDens,const std::string nFPO2Stat,
                     const std::string nFVegfStat);
void toyModel(double *x, double *y);

#endif
