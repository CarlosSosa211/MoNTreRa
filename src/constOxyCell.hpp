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
#define CONSTOXYCELL_NUM_IN_B 3
#define CONSTOXYCELL_NUM_IN_I 0
#define CONSTOXYCELL_NUM_IN_D 1

#define IN_OXY_PO2  m_inD[0] //pO2 (mmHg)

//State variables
#define CONSTOXYCELL_NUM_ST_B 7
#define CONSTOXYCELL_NUM_ST_I 0
#define CONSTOXYCELL_NUM_ST_D 2

//Outputs
#define CONSTOXYCELL_NUM_OUT_B 0
#define CONSTOXYCELL_NUM_OUT_I 0
#define CONSTOXYCELL_NUM_OUT_D 2

//Internal parameters
#define CONSTOXYCELL_NUM_PAR_B 0
#define CONSTOXYCELL_NUM_PAR_I 0
#define CONSTOXYCELL_NUM_PAR_D 4

class ConstOxyCell : public AbsOxyCell{
public :
    ConstOxyCell(const double hypThres, Model *const parent);
    virtual ~ConstOxyCell();
    virtual int updateModel(const double currentTime, const double DT);
    void setInPO2(const double input);
};

#endif
