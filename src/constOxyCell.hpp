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

//Inputs
#define IN_OXYPO2      m_in->at(3)
#define IN_OXYVEGF     m_in->at(4)

class ConstOxyCell : public AbsOxyCell{
public :
    ConstOxyCell();
    ConstOxyCell(const double pO2NotVes, const double pO2NormVes,
                 const double pO2TumVes, const double hypThres,
                 Model *const parent);
    ConstOxyCell(const double hypThres, Model *const parent);
    virtual ~ConstOxyCell();
    virtual int updateModel(const double currentTime, const double DT);
    void setInPO2(const double input);
};

#endif
