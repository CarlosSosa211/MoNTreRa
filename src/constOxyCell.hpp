/**
 * @file constOxyCell.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.07.17
 */

#ifndef DEF_CONSTOXYCELL
#define DEF_CONSTOXYCELL

#include "absOxyCell.hpp"

class ConstOxyCell : public AbsOxyCell{
public :
    ConstOxyCell();
    ConstOxyCell(const double pO2NotVes, const double pO2NormVes,
                 const double pO2TumVes, const double hypThres,
                 Model *const parent);
    virtual ~ConstOxyCell();
    virtual int updateModel(const double currentTime, const double DT);
};

#endif
