/**
 * @file constOxyTissue.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.07.17
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>

#include "constOxyTissue.hpp"

using namespace std;

/*------------------------------------------------------------------------------
 * Constructor of the class ConstOxyTissue.
 *
 * Inputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - oxy: oxygenation scenario (0, no oxygenation; 2, space dependent;
 *  3, time dependent, 4 constant),
 *  - hypThres: pO2 hypoxia threshold (mmHg).
------------------------------------------------------------------------------*/

ConstOxyTissue::ConstOxyTissue(const int nrow, const int ncol, const int nlayer,
                               const vector<bool> &inVes, const int oxy,
                               const double pO2NormVes, const double pO2TumVes,
                               const double hypThres, const double pO2NotVes) :
    AbsOxyTissue(CONSTOXYTISSUE_NUM_IN_B, CONSTOXYTISSUE_NUM_IN_I,
                 CONSTOXYTISSUE_NUM_IN_D, CONSTOXYTISSUE_NUM_ST_B,
                 CONSTOXYTISSUE_NUM_ST_I, CONSTOXYTISSUE_NUM_ST_D,
                 CONSTOXYTISSUE_NUM_OUT_B, CONSTOXYTISSUE_NUM_OUT_I,
                 CONSTOXYTISSUE_NUM_OUT_D, CONSTOXYTISSUE_NUM_PAR_B,
                 CONSTOXYTISSUE_NUM_PAR_I, CONSTOXYTISSUE_NUM_PAR_D, nrow *
                 ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new ConstOxyCell(pO2NormVes, pO2TumVes, hypThres, this);
        ((ConstOxyCell *)m_comp->at(k))->setInNormVes(inVes.at(k));
    }

    PAR_OXYTISSUE_OXY = oxy;
    PAR_PO2_NOT_VES   = pO2NotVes;
}


/*------------------------------------------------------------------------------
 * Destructor of the class ConstOxyTissue.
------------------------------------------------------------------------------*/

ConstOxyTissue::~ConstOxyTissue(){
    for(int k(0); k < m_numComp; k++){
        delete m_comp->at(k);
    }
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model initModel method.
------------------------------------------------------------------------------*/

int ConstOxyTissue::initModel(){
    for (int k(0); k < m_numComp; k++){
        (m_comp->at(k))->initModel();
    }
    return 1;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model calcModelOut method.
------------------------------------------------------------------------------*/

int ConstOxyTissue::calcModelOut(){
    OUT_HYP_DENS = double(getNumHyp()) / double(m_numComp) * 100.0;

    vector<double> pO2;
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->calcModelOut();
        pO2.push_back(m_comp->at(k)->getOutD()[0]);
    }
    int n (pO2.size() / 2);
    nth_element(pO2.begin(), pO2.begin() + n, pO2.end());
    OUT_PO2_MED = pO2.at(n);

    OUT_PO2_MEAN = accumulate(pO2.begin(), pO2.end(), 0.0) / pO2.size();

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model updateModel method.
 *
 * Inputs:
 *  - currentTime: simulation current time (ms),
 *  - DT: simulation timestep (ms).
------------------------------------------------------------------------------*/

int ConstOxyTissue::updateModel(double currentTime, const double DT){
    double _numComp100(1.0 / double(m_numComp) * 100.0);
    ST_OXYTISSUE_VES_DENS      = double(getNumVes()) * _numComp100;
    ST_OXYTISSUE_NORM_VES_DENS = double(getNumNormVes()) * _numComp100;
    ST_OXYTISSUE_TUM_VES_DENS  = double(getNumTumVes()) * _numComp100;

    ST_OXYTISSUE_DEAD_DENS  = double(getNumDead()) * _numComp100;

    switch(PAR_OXYTISSUE_OXY){
    case 0:{
        for(int k(0); k < m_numComp; k++){
            ((ConstOxyCell *)m_comp->at(k))->setInPO2(PAR_PO2_NOT_VES);
        }
        break;
    }
    case 3:{
        for(int k(0); k < m_numComp; k++){
            ((ConstOxyCell *)m_comp->at(k))->setInPO2(0.1 *
                                                      ST_OXYTISSUE_DEAD_DENS);
        }
        break;
    }
    case 4:{
        for(int k(0); k < m_numComp; k++){
            ((ConstOxyCell *)m_comp->at(k))->setInPO2(PAR_PO2_NOT_VES);
        }
        break;
    }
    }

    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->updateModel();
    }
    return 1;
}
