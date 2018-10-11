/**
 * @file constOxyCell.cpp
 * @brief
 * @author Nicolas Ciferri
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.07.17
 */

#include <iostream>
#include "constOxyCell.hpp"

using namespace std;

ConstOxyCell::ConstOxyCell() : AbsOxyCell(){
}


ConstOxyCell::ConstOxyCell(const double pO2, const double pO2NormVes,
                           const double pO2TumVes, const double hypThres,
                           Model *const parent) : AbsOxyCell(){
    PAR_OXYPO2       = pO2;
    PAR_PO2_NORM_VES = pO2NormVes;
    PAR_PO2_TUM_VES  = pO2TumVes;
    PAR_HYP_THRES    = hypThres;
    m_parent = parent;
}


ConstOxyCell::~ConstOxyCell(){
}


int ConstOxyCell::updateModel(const double currentTime,
                              const double DT){
    ST_OXYDEAD     = IN_OXYDEAD;
    ST_OXYNORM_VES = IN_OXYNORM_VES;
    ST_OXYTUM_VES  = IN_OXYTUM_VES;

    if(ST_OXYNORM_VES){
        ST_PO2 = PAR_PO2_NORM_VES;
    }
    else if(ST_OXYTUM_VES){
        ST_PO2 = PAR_PO2_TUM_VES;
    }
    else if(ST_OXYDEAD){
        ST_PO2 = 0.0;
    }
    else{
        ST_PO2 = PAR_OXYPO2;
    }
    return 0;
}




