#ifndef DEF_SENSANR
#define DEF_SENSANR

#include <fstream>
#include <iostream>

void evalR(const int nMethod, const int nModel,
           const std::string nFInTissueDim = "", const std::string nFInTum = "",
           const std::string nFInVes = "", const std::string nFInPO2 = "");

#endif
