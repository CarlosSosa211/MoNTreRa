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

ConstOxyCell::ConstOxyCell(const double hypThres, Model *const parent) :
    AbsOxyCell(){

    m_in->resize(5);

    PAR_HYP_THRES = hypThres;
    m_parent = parent;
}


ConstOxyCell::~ConstOxyCell(){
}


/*int ConstOxyCell::initModel(){
    ST_OXYDEAD     = IN_OXYDEAD;
    ST_OXYNORM_VES = IN_OXYNORM_VES;
    ST_OXYTUM_VES  = IN_OXYTUM_VES;

    ST_OXYPO2  = IN_OXYPO2;
    ST_OXYVEGF = IN_OXYVEGF;
    ST_HYP = 0.0;

    ST_OXYSTABLE_CELL  = 0.0;
    ST_VEGFSTABLE_CELL = 0.0;

    return 0;
}*/


int ConstOxyCell::updateModel(const double currentTime, const double DT){
    ST_OXYDEAD     = IN_OXYDEAD;
    ST_OXYNORM_VES = IN_OXYNORM_VES;
    ST_OXYTUM_VES  = IN_OXYTUM_VES;

    ST_OXYPO2  = IN_OXYPO2;
    ST_OXYVEGF = IN_OXYVEGF;

    ST_HYP = !ST_OXYDEAD * (ST_OXYPO2 < PAR_HYP_THRES);

    ST_OXYSTABLE_CELL  = 1.0;
    ST_VEGFSTABLE_CELL = 1.0;

    return 0;
}


void ConstOxyCell::setInPO2(const double input){
    IN_OXYPO2 = input;
}




