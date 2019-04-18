/**
 * @file constOxyTissue.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.07.17
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>

#include "constOxyTissue.hpp"

using namespace std;


ConstOxyTissue::ConstOxyTissue(const int nrow, const int ncol, const int nlayer,
                               const vector<bool> &inVes,
                               const vector<double> &inPO2,
                               const double hypThres) : Model(0, 0, 5, 0, nrow *
                                                              ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new ConstOxyCell(hypThres, this);
        m_numOut += (m_comp->at(k))->getNumOut();
        ((ConstOxyCell *)m_comp->at(k))->setInNormVes(inVes.at(k));
        ((ConstOxyCell *)m_comp->at(k))->setInPO2(inPO2.at(k));
    }
}


ConstOxyTissue::~ConstOxyTissue(){
    for(int k(0); k < m_numComp; k++){
        delete m_comp->at(k);
    }
}


int ConstOxyTissue::initModel(){
    for (int k(0); k < m_numComp; k++){
        (m_comp->at(k))->initModel();
    }
    return 1;
}


int ConstOxyTissue::calcModelOut(){
    OUT_HYP_DENS = double(getNumHyp()) / double(m_numComp) * 100.0;

    vector<double> pO2;
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->calcModelOut();
        pO2.push_back(m_comp->at(k)->getOut()->at(0));
    }
    int n (pO2.size() / 2);
    nth_element(pO2.begin(), pO2.begin() + n, pO2.end());
    OUT_PO2_MED = pO2.at(n);

    OUT_PO2_MEAN = accumulate(pO2.begin(), pO2.end(), 0.0) / pO2.size();

    return 0;
}


int ConstOxyTissue::updateModel(double currentTime, const double DT){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->updateModel();
    }
    return 1;
}


int ConstOxyTissue::terminateModel(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->terminateModel();
    }
    return 0;
}


int ConstOxyTissue::getNumHyp() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((ConstOxyCell *)m_comp->at(k))->getHyp()){
            count++;
        }
    }
    return count;
}
