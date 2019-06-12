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
 *  - pO2NormVes: fixed pO2 value for pre-existing endothelial cells (mmHg),
 *  - pO2TumVes: fixed pO2 value for neo-created endothelial cells (mmHg),
 *  - hypThres: pO2 hypoxia threshold (mmHg),
 *  - parent: pointer to the parent of the ConstOxyCell, a ConstOxyTissue.
------------------------------------------------------------------------------*/

ConstOxyCell::ConstOxyCell(const double pO2NormVes, const double pO2TumVes,
                           const double hypThres, Model *const parent) :
    AbsOxyCell(pO2NormVes, pO2TumVes, hypThres, CONSTOXYCELL_NUM_IN_B,
               CONSTOXYCELL_NUM_IN_I, CONSTOXYCELL_NUM_IN_D,
               CONSTOXYCELL_NUM_ST_B, CONSTOXYCELL_NUM_ST_I,
               CONSTOXYCELL_NUM_ST_D, CONSTOXYCELL_NUM_OUT_B,
               CONSTOXYCELL_NUM_OUT_I, CONSTOXYCELL_NUM_OUT_D,
               CONSTOXYCELL_NUM_PAR_B, CONSTOXYCELL_NUM_PAR_I,
               CONSTOXYCELL_NUM_PAR_D){
    IN_OXY_PO2 = 0.0;
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
    else{
        ST_OXYCELL_PO2 = IN_OXY_PO2;
    }

    ST_HYP = !ST_OXYCELL_DEAD && (ST_OXYCELL_PO2 < PAR_HYP_THRES);

    ST_OXYCELL_OXY_STABLE  = true;
    ST_OXYCELL_VEGF_STABLE = true;

    return 0;
}

/*------------------------------------------------------------------------------
 * This function sets the pO2 input.
 *
 * Inputs:
 *  - input: pO2 input.
------------------------------------------------------------------------------*/

void ConstOxyCell::setInPO2(const double input){
    IN_OXY_PO2 = input;
}
