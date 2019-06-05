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

/*------------------------------------------------------------------------------
 * Constructor of the class Cell.
------------------------------------------------------------------------------*/

AbsOxyCell::AbsOxyCell(const int numInB, const int numInI, const int numInD,
                       const int numStB, const int numStI, const int numStD,
                       const int numOutB, const int numOutI, const int numOutD,
                       const int numParB, const int numParI,
                       const int numParD) :
    Model(numInB, numInI, numInD, numStB, numStI, numStD, numOutB, numOutI,
          numOutD, numParB, numParI, numParD, 0){
}


/*------------------------------------------------------------------------------
 * Destructor of the class Cell.
------------------------------------------------------------------------------*/

AbsOxyCell::~AbsOxyCell(){
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model calcModelOut method.
------------------------------------------------------------------------------*/

int AbsOxyCell::calcModelOut(){
    OUT_PO2  = ST_OXYCELL_PO2;
    OUT_VEGF = ST_OXYCELL_VEGF;

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model initModel method.
------------------------------------------------------------------------------*/

int AbsOxyCell::initModel(){
    ST_OXYCELL_DEAD     = IN_OXYCELL_DEAD;
    ST_OXYCELL_NORM_VES = IN_OXYCELL_NORM_VES;
    ST_OXYCELL_TUM_VES  = IN_OXYCELL_TUM_VES;
    ST_OXYCELL_VES      = ST_OXYCELL_NORM_VES || ST_OXYCELL_TUM_VES;

    if(ST_OXYCELL_NORM_VES){
        ST_OXYCELL_PO2 = PAR_PO2_NORM_VES;
    }
    else if(ST_OXYCELL_TUM_VES){
        ST_OXYCELL_PO2 = PAR_PO2_TUM_VES;
    }
    else if(ST_OXYCELL_DEAD){
        ST_OXYCELL_PO2 = 0.0;
    }
    else{
        ST_OXYCELL_PO2 = 10.0;
    }

    ST_HYP = false;
    ST_OXYCELL_VEGF = 0.0;

    ST_OXYCELL_OXY_STABLE  = false;
    ST_OXYCELL_VEGF_STABLE = false;

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model terminalModel method.
------------------------------------------------------------------------------*/

int AbsOxyCell::terminateModel(){
    return 0;
}


/*------------------------------------------------------------------------------
 * This function gets the healthy cell state.
------------------------------------------------------------------------------*/

bool AbsOxyCell::getDead() const{
    return ST_OXYCELL_DEAD;
}


/*------------------------------------------------------------------------------
 * This function gets the hypoxic state.
------------------------------------------------------------------------------*/

bool AbsOxyCell::getHyp() const{
    return ST_HYP;
}


/*------------------------------------------------------------------------------
 * This function gets the pre-existing endothelial cell state.
------------------------------------------------------------------------------*/

bool AbsOxyCell::getNormVes() const{
    return ST_OXYCELL_NORM_VES;
}


/*------------------------------------------------------------------------------
 * This function gets the pO2 output.
------------------------------------------------------------------------------*/
double AbsOxyCell::getOutPO2() const{
    return OUT_PO2;
}


/*------------------------------------------------------------------------------
 * This function gets the VEGF concentration output.
------------------------------------------------------------------------------*/

double AbsOxyCell::getOutVEGF() const{
    return OUT_VEGF;
}


/*------------------------------------------------------------------------------
 * This function gets the pO2 stable state.
------------------------------------------------------------------------------*/

bool AbsOxyCell::getOxyStable() const{
    return ST_OXYCELL_OXY_STABLE;
}


/*------------------------------------------------------------------------------
 * This function gets the pO2 state.
------------------------------------------------------------------------------*/

double AbsOxyCell::getPO2() const{
    return ST_OXYCELL_PO2;
}


/*------------------------------------------------------------------------------
 * This function gets the neo-created endothelial cell state.
------------------------------------------------------------------------------*/

bool AbsOxyCell::getTumVes() const{
    return ST_OXYCELL_TUM_VES;
}


/*------------------------------------------------------------------------------
 * This function gets the VEGF concentration state.
------------------------------------------------------------------------------*/

double AbsOxyCell::getVEGF() const{
    return ST_OXYCELL_VEGF;
}


/*------------------------------------------------------------------------------
 * This function gets the VEGF concentration stable state.
------------------------------------------------------------------------------*/

bool AbsOxyCell::getVegfStable() const{
    return ST_OXYCELL_VEGF_STABLE;
}


/*------------------------------------------------------------------------------
 * This function gets the endothelial cell state.
------------------------------------------------------------------------------*/

bool AbsOxyCell::getVes() const{
    return ST_OXYCELL_VES;
}


/*------------------------------------------------------------------------------
 * This function sets the dead cell input.
 *
 * Inputs:
 *  - input: dead cell input.
------------------------------------------------------------------------------*/

void AbsOxyCell::setInDead(const bool input){
    IN_OXYCELL_DEAD = input;
}


/*------------------------------------------------------------------------------
 * This function sets the pre-existing endothelial cell input.
 *
 * Inputs:
 *  - input: pre-existing endothelial cell input.
------------------------------------------------------------------------------*/

void AbsOxyCell::setInNormVes(const bool input){
    IN_OXYCELL_NORM_VES = input;
}


/*------------------------------------------------------------------------------
 * This function sets the neo-created endothelial cell input.
 *
 * Inputs:
 *  - input: neo-created endothelial cell input.
------------------------------------------------------------------------------*/

void AbsOxyCell::setInTumVes(const bool input){
    IN_OXYCELL_TUM_VES = input;
}
