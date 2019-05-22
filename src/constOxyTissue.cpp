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
 *  - oxy: integer indicating the oxygenation scenario (0, no oxygenation;
 *  2, space dependent; 3, time dependent, 4 constant),
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - inPO2: vector containing the initial pO2 values,
 *  - hypThres: pO2 hypoxia threshold (mmHg).
------------------------------------------------------------------------------*/

ConstOxyTissue::ConstOxyTissue(const int nrow, const int ncol, const int nlayer,
                               const vector<bool> &inVes,
                               const vector<double> &inPO2,
                               const double hypThres) :
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
        m_comp->at(k) = new ConstOxyCell(hypThres, this);
        ((ConstOxyCell *)m_comp->at(k))->setInNormVes(inVes.at(k));
        ((ConstOxyCell *)m_comp->at(k))->setInPO2(inPO2.at(k));
    }
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
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->updateModel();
    }
    return 1;
}
