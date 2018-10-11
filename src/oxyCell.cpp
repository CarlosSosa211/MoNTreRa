/**
 * @file oxyCell.cpp
 * @brief
 * @author Nicolas Ciferri
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include "oxyCell.hpp"

using namespace std;

OxyCell::OxyCell() : AbsOxyCell(){
}


OxyCell::OxyCell(const double Vmax, const double Km,
                 const double pO2NormVes, const double pO2TumVes,
                 const double hypThres, const double VmaxVegf,
                 const double KmVegf, const double hypVegf,
                 Model *const parent) : AbsOxyCell(){
    m_in->resize(7);
    m_param->resize(9);

    PAR_VMAX         = Vmax;
    PAR_KM           = Km;
    PAR_OXYPO2       = 2.0 * hypThres;
    PAR_PO2_NORM_VES = pO2NormVes;
    PAR_PO2_TUM_VES  = pO2TumVes;
    PAR_HYP_THRES    = hypThres;
    PAR_VMAX_VEGF    = VmaxVegf;
    PAR_KM_VEGF      = KmVegf;
    PAR_HYP_VEGF     = hypVegf;
    m_parent = parent;
}


OxyCell::~OxyCell(){
}


int OxyCell::updateModel(const double currentTime,
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
    else{
        calcConsO2();
        ST_PO2 += IN_DIFF_O2 - IN_CONS_O2;
    }
    ST_HYP = !ST_OXYDEAD * (ST_PO2 < PAR_HYP_THRES);
    if(ST_HYP){
        ST_VEGF = PAR_HYP_VEGF;
    }
    else{
        calcConsVegf();
        ST_VEGF += IN_DIFF_VEGF - IN_CONS_VEGF;
    }

    return 0;
}


void OxyCell::calcConsO2(){
    if(ST_OXYDEAD){
        IN_CONS_O2 = 0.0;
    }
    else{
        IN_CONS_O2 = ST_PO2 * PAR_VMAX /
                (PAR_KM + ST_PO2);
    }
}


void OxyCell::calcConsVegf(){
    if(ST_OXYNORM_VES || ST_OXYTUM_VES){
        IN_CONS_VEGF = ST_VEGF * PAR_VMAX_VEGF /
                (PAR_KM_VEGF + ST_VEGF);
    }
    else{
        IN_CONS_VEGF = 0.0;
    }
}


void OxyCell::setInDiffO2(const double input){
    IN_DIFF_O2 = input;
}


void OxyCell::setInDiffVEGF(const double input){
    IN_DIFF_VEGF = input;
}




