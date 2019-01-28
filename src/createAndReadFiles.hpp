#include <vector>

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double sigmaTum,
                  const double vascDens, const double sigmaVasc,
                  std::vector<bool> &inTum, std::vector<bool> &inVes);

void readInFiles(const std::string nFInTissueDim, const std::string nFInTum,
                 const std::string nFInVes, int &nrow, int &ncol, int &nlayer,
                 double &cellSize, std::vector<bool> &inTum, std::vector<bool> &inVes);
