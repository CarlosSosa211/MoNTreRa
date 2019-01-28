#include <string>
#include <vector>

#include "constOxyTissue.hpp"
#include "createAndReadFiles.hpp"
#include "oxyTissue.hpp"
#include "rootSimulator.hpp"
#include "simpcell.hpp"
#include "tissue.hpp"
#include "treatment.hpp"

void tcp(const int N, const std::string nFInTissueTCP, const std::string nFParTCP,
         const std::vector<std::string> nFTreatmentTCP);
void modelTCP(const double *x, double *y, const int nrow,
              const int ncol, const int nlayer, const double cellSize,
              const std::vector<bool> &inTum, const std::vector<bool> & inVes,
              Treatment *const treatment);
