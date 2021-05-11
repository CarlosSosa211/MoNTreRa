#ifndef DEF_SENSANR
#define DEF_SENSANR

#include <fstream>
#include <iostream>

void evalPy(const int nMethod, const std::string nFX, const std::string nFY);
void evalPy(const int nMethod, const std::string nFInTissuePar,
            const std::string nFX, const std::string nFTreatment,
            const std::string nFY);

#endif
