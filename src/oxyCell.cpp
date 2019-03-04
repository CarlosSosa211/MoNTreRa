/**
 * @file oxyCell.cpp
 * @brief
 * @author Nicolas Ciferri
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include <iostream>

#include "oxyCell.hpp"

using namespace std;

OxyCell::OxyCell() : AbsOxyCell(){
}


OxyCell::OxyCell(const double ang, const double  VmaxVegf, const double KmVegf,
                 const double hypVegf, const double VmaxO2, const double KmO2,
                 const double pO2NormVes, const double pO2TumVes,
                 const double hypThres, Model *const parent) : AbsOxyCell(){
    m_in->resize(7);
    m_param->resize(10);

    PAR_OXYCELL_ANG  = ang;
    PAR_VMAX_VEGF    = VmaxVegf;
    PAR_KM_VEGF      = KmVegf;
    PAR_HYP_VEGF     = hypVegf;

    PAR_VMAX_O2      = VmaxO2;
    PAR_KM_O2        = KmO2;
    PAR_PO2_NORM_VES = pO2NormVes;
    PAR_PO2_TUM_VES  = pO2TumVes;
    PAR_HYP_THRES    = hypThres;
    m_parent = parent;

    m_edge = new vector<OxyCell *>((unsigned int)0, 0);
}


OxyCell::~OxyCell(){
}


int OxyCell::updateModel(const double currentTime, const double DT){
    ST_OXYDEAD     = IN_OXYDEAD;
    ST_OXYNORM_VES = IN_OXYNORM_VES;
    ST_OXYTUM_VES  = IN_OXYTUM_VES;

    if(ST_OXYNORM_VES){
        ST_OXYPO2 = PAR_PO2_NORM_VES;
    }
    else if(ST_OXYTUM_VES){
        ST_OXYPO2 = PAR_PO2_TUM_VES;
    }
    else{
        calcConsO2();
        ST_OXYPO2 += IN_DIFF_O2 - IN_CONS_O2;
    }

    if(fabs(ST_OXYPO2 - OUT_PO2) < 1e-2){
        ST_OXYSTABLE_CELL = 1.0;
    }
    else{
        ST_OXYSTABLE_CELL = 0.0;
    }

    ST_HYP = !ST_OXYDEAD * (ST_OXYPO2 < PAR_HYP_THRES);

    if(PAR_OXYCELL_ANG){
        if(ST_HYP){
            ST_OXYVEGF = PAR_HYP_VEGF;
        }
        else{
            calcConsVegf();
            ST_OXYVEGF += IN_DIFF_VEGF - IN_CONS_VEGF;
        }
    }

    if(fabs(ST_OXYVEGF - OUT_VEGF) < 1e-2){
        ST_VEGFSTABLE_CELL = 1.0;
    }
    else{
        ST_VEGFSTABLE_CELL = 0.0;
    }

    return 0;
}


void OxyCell::addToEdge(OxyCell *const cell){
    m_edge->push_back(cell);
}


void OxyCell::calcConsO2(){
    if(ST_OXYDEAD){
        IN_CONS_O2 = 0.0;
    }
    else{
        IN_CONS_O2 = ST_OXYPO2 * PAR_VMAX_O2 / (PAR_KM_O2 + ST_OXYPO2);
    }
}


void OxyCell::calcConsVegf(){
    if(ST_OXYNORM_VES || ST_OXYTUM_VES){
        IN_CONS_VEGF = ST_OXYVEGF * PAR_VMAX_VEGF / (PAR_KM_VEGF + ST_OXYVEGF);
    }
    else{
        IN_CONS_VEGF = 0.0;
    }
}


vector<OxyCell *> *OxyCell::getEdge() const{
    return m_edge;
}


void OxyCell::setInDiffO2(const double input){
    IN_DIFF_O2 = input;
}


void OxyCell::setInDiffVEGF(const double input){
    IN_DIFF_VEGF = input;
}




