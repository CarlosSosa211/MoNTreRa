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

/*------------------------------------------------------------------------------
 * Constructor of the class ConstOxyCell.
 *
 * Inputs:
 *  - hypThres: pO2 hypoxia threshold (mmHg),
 *  - parent: pointer to the parent of the ConstOxyCell, a ConstOxyTissue.
------------------------------------------------------------------------------*/

ConstOxyCell::ConstOxyCell(const double hypThres, Model *const parent) :
    AbsOxyCell(CONSTOXYCELL_NUM_IN_B, CONSTOXYCELL_NUM_IN_I,
               CONSTOXYCELL_NUM_IN_D, CONSTOXYCELL_NUM_ST_B,
               CONSTOXYCELL_NUM_ST_I, CONSTOXYCELL_NUM_ST_D,
               CONSTOXYCELL_NUM_OUT_B, CONSTOXYCELL_NUM_OUT_I,
               CONSTOXYCELL_NUM_OUT_D, CONSTOXYCELL_NUM_PAR_B,
               CONSTOXYCELL_NUM_PAR_I, CONSTOXYCELL_NUM_PAR_D){
    PAR_HYP_THRES = hypThres;
    m_parent = parent;
}


/*------------------------------------------------------------------------------
 * Destructor of the class ConstOxyCell.
------------------------------------------------------------------------------*/

ConstOxyCell::~ConstOxyCell(){
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model updateModel method.
 *
 * Inputs:
 *  - currentTime: simulation current time (ms),
 *  - DT: simulation timestep (ms).
------------------------------------------------------------------------------*/

int ConstOxyCell::updateModel(const double currentTime, const double DT){
    ST_OXY_DEAD     = IN_OXY_DEAD;
    ST_OXY_NORM_VES = IN_OXY_NORM_VES;
    ST_OXY_TUM_VES  = IN_OXY_TUM_VES;

    ST_OXY_PO2  = IN_OXY_PO2;

    ST_HYP = !ST_OXY_DEAD && (ST_OXY_PO2 < PAR_HYP_THRES);

    ST_OXY_STABLE_CELL  = true;
    ST_VEGF_STABLE_CELL = true;

    return 0;
}


/*------------------------------------------------------------------------------
 * This function sets the pO2 input.
 *
 * Inputs:
 *  - input: pO2 input (mmHg).
------------------------------------------------------------------------------*/

void ConstOxyCell::setInPO2(const double input){
    IN_OXY_PO2 = input;
}
