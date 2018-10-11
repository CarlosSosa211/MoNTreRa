/**
 * @file absOxyCell.cpp
 * @brief
 * @author Nicolas Ciferri
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.07.17
 */

#include <iostream>
#include "absOxyCell.hpp"

using namespace std;

AbsOxyCell::AbsOxyCell() : Model(3, 6, 2, 5, 0){
}


AbsOxyCell::~AbsOxyCell(){
}


int AbsOxyCell::calcModelOut(){
    OUT_PO2  = ST_PO2;
    OUT_VEGF = ST_VEGF;

    return 0;
}


int AbsOxyCell::initModel(){
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
        ST_HYP = ST_PO2 < PAR_HYP_THRES;
    }

    if(ST_HYP){
        ST_VEGF = PAR_HYP_VEGF;
    }
    else{
        ST_VEGF = 0.0;
    }
    return 0;
}


int AbsOxyCell::terminateModel(){
    return 0;
}


bool AbsOxyCell::getHyp() const{
    return ST_HYP;
}


double AbsOxyCell::getOutPO2() const{
    return OUT_PO2;
}


double AbsOxyCell::getOutVEGF() const{
    return OUT_VEGF;
}


double AbsOxyCell::getPO2() const{
    return ST_PO2;
}


double AbsOxyCell::getVEGF() const{
    return ST_VEGF;
}


void AbsOxyCell::setInDead(const double input){
    IN_OXYDEAD = input;
}


void AbsOxyCell::setInNormVes(const double input){
    IN_OXYNORM_VES = input;
}


void AbsOxyCell::setInTumVes(const double input){
    IN_OXYTUM_VES = input;
}



