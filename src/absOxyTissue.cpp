/**
 * @file AbsOxyTissue.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.22.19
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>

#include "absOxyTissue.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * Constructor of the class AbsOxyTissue.
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

AbsOxyTissue::AbsOxyTissue(const int numInB, const int numInI, const int numInD,
                           const int numStB, const int numStI, const int numStD,
                           const int numOutB, const int numOutI,
                           const int numOutD, const int numParB,
                           const int numParI, const int numParD,
                           const int numComp) :
    Model(numInB, numInI, numInD, numStB, numStI, numStD, numOutB, numOutI,
          numOutD, numParB, numParI, numParD, numComp){
}


/*------------------------------------------------------------------------------
 * Destructor of the class AbsOxyTissue.
------------------------------------------------------------------------------*/

AbsOxyTissue::~AbsOxyTissue(){
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model initModel method.
------------------------------------------------------------------------------*/

int AbsOxyTissue::initModel(){
    for (int k(0); k < m_numComp; k++){
        (m_comp->at(k))->initModel();
    }
    return 1;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model terminalModel method.
------------------------------------------------------------------------------*/

int AbsOxyTissue::terminateModel(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->terminateModel();
    }
    return 0;
}


/*------------------------------------------------------------------------------
 * This function counts the hypoxic cells of the tissue.
 *
 * Outputs:
 *  - count: number of hypoxic cells of the tissue
------------------------------------------------------------------------------*/

int AbsOxyTissue::getNumHyp() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((AbsOxyCell *)m_comp->at(k))->getHyp()){
            count++;
        }
    }
    return count;
}
