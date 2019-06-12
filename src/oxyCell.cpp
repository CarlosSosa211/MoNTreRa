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

/*------------------------------------------------------------------------------
 * Constructor of the class OxyCell.
 *
 * Inputs:
 *  - ang: angiogenesis,
 *  - VmaxVegf: maximum VEGF consumption ratio (mol/um^3ms),
 *  - KmVegf: VEGF Michaelis constant (mol/um^3),
 *  - hypVegf: VEGF fixed value for hypoxic cells (mol/um^3),
 *  - VmaxO2: maximum pO2 consumption ratio (mmHg/ms),
 *  - KmO2: pO2 Michaelis constant (mmHg),
 *  - pO2NormVes: fixed pO2 value for pre-existing endothelial cells (mmHg),
 *  - pO2TumVes: fixed pO2 value for neo-created endothelial cells (mmHg),
 *  - hypThres: pO2 hypoxia threshold (mmHg),
 *  - parent: pointer to the parent of the Oxyell, an OxyTissue.
------------------------------------------------------------------------------*/

OxyCell::OxyCell(const bool ang, const double  VmaxVegf, const double KmVegf,
                 const double hypVegf, const double VmaxO2, const double KmO2,
                 const double pO2NormVes, const double pO2TumVes,
                 const double hypThres, Model *const parent) :
    AbsOxyCell(pO2NormVes, pO2TumVes, hypThres, OXYCELL_NUM_IN_B,
               OXYCELL_NUM_IN_I, OXYCELL_NUM_IN_D, OXYCELL_NUM_ST_B,
               OXYCELL_NUM_ST_I, OXYCELL_NUM_ST_D, OXYCELL_NUM_OUT_B,
               OXYCELL_NUM_OUT_I, OXYCELL_NUM_OUT_D, OXYCELL_NUM_PAR_B,
               OXYCELL_NUM_PAR_I, OXYCELL_NUM_PAR_D){
    IN_DIFF_O2   = 0.0;
    IN_CONS_O2   = 0.0;
    IN_DIFF_VEGF = 0.0;
    IN_CONS_VEGF = 0.0;

    PAR_OXYCELL_ANG  = ang;
    PAR_VMAX_VEGF    = VmaxVegf;
    PAR_KM_VEGF      = KmVegf;
    PAR_HYP_VEGF     = hypVegf;

    PAR_VMAX_O2      = VmaxO2;
    PAR_KM_O2        = KmO2;
    m_parent = parent;

    m_edge = new vector<OxyCell *>((unsigned int)0, 0);
}


/*------------------------------------------------------------------------------
 * Destructor of the class OxyCell.
------------------------------------------------------------------------------*/

OxyCell::~OxyCell(){
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model updateModel method.
 *
 * Inputs:
 *  - currentTime: simulation current time (ms),
 *  - DT: simulation timestep (ms).
------------------------------------------------------------------------------*/

int OxyCell::updateModel(const double currentTime, const double DT){
    ST_OXYCELL_DEAD     = IN_OXYCELL_DEAD;
    ST_OXYCELL_NORM_VES = IN_OXYCELL_NORM_VES;
    ST_OXYCELL_TUM_VES  = IN_OXYCELL_TUM_VES;

    if(ST_OXYCELL_NORM_VES){
        ST_OXYCELL_PO2 = PAR_PO2_NORM_VES;
    }
    else if(ST_OXYCELL_TUM_VES){
        ST_OXYCELL_PO2 = PAR_PO2_TUM_VES;
    }
    else{
        calcConsO2();
        ST_OXYCELL_PO2 += IN_DIFF_O2 - IN_CONS_O2;
    }

    if(fabs(ST_OXYCELL_PO2 - OUT_PO2) < 1e-2){
        ST_OXYCELL_OXY_STABLE = true;
    }
    else{
        ST_OXYCELL_OXY_STABLE = false;
    }

    ST_HYP = !ST_OXYCELL_DEAD && (ST_OXYCELL_PO2 < PAR_HYP_THRES);

    if(PAR_OXYCELL_ANG){
        if(ST_HYP){
            ST_OXYCELL_VEGF = PAR_HYP_VEGF;
        }
        else{
            calcConsVegf();
            ST_OXYCELL_VEGF += IN_DIFF_VEGF - IN_CONS_VEGF;
        }
    }

    if(fabs(ST_OXYCELL_VEGF - OUT_VEGF) < 1e-2){
        ST_OXYCELL_VEGF_STABLE = true;
    }
    else{
        ST_OXYCELL_VEGF_STABLE = false;
    }

    return 0;
}


/*------------------------------------------------------------------------------
 * This function adds a cell to the edge of the current one.
 *
 * Inputs:
 *  - cell: pointer to the cell to be added to the edge of the current one.
------------------------------------------------------------------------------*/

void OxyCell::addToEdge(OxyCell *const cell){
    m_edge->push_back(cell);
}


/*------------------------------------------------------------------------------
 * This function calculates the pO2 consumption of the current cell.
------------------------------------------------------------------------------*/

void OxyCell::calcConsO2(){
    if(ST_OXYCELL_DEAD){
        IN_CONS_O2 = 0.0;
    }
    else{
        IN_CONS_O2 = ST_OXYCELL_PO2 * PAR_VMAX_O2 / (PAR_KM_O2 +
                                                     ST_OXYCELL_PO2);
    }
}


/*------------------------------------------------------------------------------
 * This function calculates the VEGF consumption of the current cell.
------------------------------------------------------------------------------*/

void OxyCell::calcConsVegf(){
    if(ST_OXYCELL_NORM_VES || ST_OXYCELL_TUM_VES){
        IN_CONS_VEGF = ST_OXYCELL_VEGF * PAR_VMAX_VEGF / (PAR_KM_VEGF +
                                                          ST_OXYCELL_VEGF);
    }
    else{
        IN_CONS_VEGF = 0.0;
    }
}


/*------------------------------------------------------------------------------
 * This function gets the edge of the current cell.
------------------------------------------------------------------------------*/

vector<OxyCell *> *OxyCell::getEdge() const{
    return m_edge;
}


/*------------------------------------------------------------------------------
 * This function sets the diffused pO2 input.
 *
 * Inputs:
 *  - input: diffused pO2 input (mmHg).
------------------------------------------------------------------------------*/

void OxyCell::setInDiffO2(const double input){
    IN_DIFF_O2 = input;
}


/*------------------------------------------------------------------------------
 * This function sets the diffused VEGF input.
 *
 * Inputs:
 *  - input: diffused VEGF input (mol/um^3).
------------------------------------------------------------------------------*/

void OxyCell::setInDiffVEGF(const double input){
    IN_DIFF_VEGF = input;
}
