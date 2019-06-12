/**
 * @file cell.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "cell.hpp"

using namespace std;

/*------------------------------------------------------------------------------
 * Constructor of the class Cell.
 *
 * Inputs:
 *  - i: row number of the cell,
 *  - j: column number of the cell,
 *  - l: layer number of the cell,
 *  - tumGrowth: tumour growth,
 *  - doubTime: duration of the cycle of tumour cells (h),
 *  - cycDur: vector containing the duration fractions of every phase of the
 *  cell cycle,
 *  - res: healthy cell division,
 *  - fibDoubTime: duration of the cycle of healthy cells (h),
 *  - ang: angiogenesis,
 *  - angTime: duration of the cycle of endothelial cells (h),
 *  - vegfThres: VEGF threshold for provoking endothelial cell division
 *  (mol/um^3),
 *  - alpha: vector containing the alpha values for every cell type and phase
 *  (Gy^-1),
 *  - beta: vector containing the beta values for every cell type and phase
 *  (Gy^-2),
 *  - doseThres: dose threshold for provoking instantaneous death by apoptosis
 *  (Gy),
 *  - arrestTime: radiation-induced arrest time (h),
 *  - oxy: oxygenation scenario (1, space-and-time dependent),
 *  - hypNecThres: pO2 hypoxic necrosis thresthold (mmHg),
 *  - parent: pointer to the parent of the Cell, an Tissue.
------------------------------------------------------------------------------*/

Cell::Cell(const int i, const int j, const int l, const bool tumGrowth,
           const double doubTime, vector <double> cycDur, const bool res,
           const double fibDoubTime, const bool ang, const double angTime,
           const double vegfThres, vector<double> alpha, vector<double> beta,
           const double doseThres, const double arrestTime, const int oxy,
           const double hypNecThres, Model *const parent) :
    Model(CELL_NUM_IN_B, CELL_NUM_IN_I, CELL_NUM_IN_D, CELL_NUM_ST_B,
          CELL_NUM_ST_I, CELL_NUM_ST_D, CELL_NUM_OUT_B, CELL_NUM_OUT_I,
          CELL_NUM_OUT_D, CELL_NUM_PAR_B, CELL_NUM_PAR_I, CELL_NUM_PAR_D, 0){
    m_i = i;
    m_j = j;
    m_l = l;

    IN_FIB      = false;
    IN_TUM      = false;
    IN_NORM_VES = false;
    IN_TUM_VES  = false;
    IN_HYP_NEC  = false;
    IN_MIT_CAT  = false;
    IN_APOP     = false;

    IN_TIMER = 0.0;
    IN_PO2   = 0.0;
    IN_VEGF  = 0.0;

    PAR_TUM_GROWTH = tumGrowth;
    PAR_DOUB_TIME  = doubTime;
    PAR_LIM_G1S    = cycDur.at(0) * PAR_DOUB_TIME;
    PAR_LIM_SG2    = PAR_LIM_G1S + cycDur.at(1) * PAR_DOUB_TIME;
    PAR_LIM_G2M    = PAR_LIM_SG2 + cycDur.at(2) * PAR_DOUB_TIME;

    PAR_RES           = res;
    PAR_FIB_DOUB_TIME = fibDoubTime;

    PAR_ANG        = ang;
    PAR_ANG_TIME   = angTime;
    PAR_VEGF_THRES = vegfThres;

    PAR_ALPHA_FIB      = alpha.at(0);
    PAR_ALPHA_G1       = alpha.at(1);
    PAR_ALPHA_S        = alpha.at(2);
    PAR_ALPHA_G2       = alpha.at(3);
    PAR_ALPHA_M        = alpha.at(4);
    PAR_ALPHA_G0       = alpha.at(5);
    PAR_ALPHA_NORM_VES = alpha.at(6);
    PAR_ALPHA_TUM_VES  = alpha.at(7);
    PAR_BETA_FIB       = beta.at(0);
    PAR_BETA_G1        = beta.at(1);
    PAR_BETA_S         = beta.at(2);
    PAR_BETA_G2        = beta.at(3);
    PAR_BETA_M         = beta.at(4);
    PAR_BETA_G0        = beta.at(5);
    PAR_BETA_NORM_VES  = beta.at(6);
    PAR_BETA_TUM_VES   = beta.at(7);
    PAR_DOSE_THRES     = doseThres;
    PAR_ARREST_TIME    = arrestTime;

    PAR_OXY = oxy;
    PAR_M = 3.0;
    PAR_K = 3.0;
    PAR_HYP_NEC_THRES = hypNecThres;

    m_parent = parent;
    m_treatment = ((Tissue *)m_parent)->getTreatment();
    m_edge = new vector<Cell *>((unsigned int)0, 0);
}


/*------------------------------------------------------------------------------
 * Destructor of the class Cell.
------------------------------------------------------------------------------*/

Cell::~Cell(){
    delete m_edge;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model calcModelOut method.
------------------------------------------------------------------------------*/

int Cell::calcModelOut(){
    OUT_STATE = ST_FIB + 2 * (ST_TUM && !ST_DAM) + 3 * (ST_TUM && ST_DAM) +
            4 * ST_NORM_VES + 5 * ST_TUM_VES + 6 * ST_HYP_NEC + 7 * ST_MIT_CAT +
            8 * ST_APOP;
    OUT_CYCLE = ST_G1 + 2 * ST_S + 3 * ST_G2 + 4 * ST_M + 5 * ST_G0;

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model initModel method.
------------------------------------------------------------------------------*/

int Cell::initModel(){
    ST_FIB      = IN_FIB;
    ST_TUM      = IN_TUM;
    ST_NORM_VES = IN_NORM_VES;
    ST_TUM_VES  = IN_TUM_VES;
    ST_VES      = ST_NORM_VES || ST_TUM_VES;
    ST_HYP_NEC  = IN_HYP_NEC;
    ST_MIT_CAT  = IN_MIT_CAT;
    ST_APOP     = IN_APOP;
    ST_DEAD     = ST_HYP_NEC || ST_MIT_CAT || ST_APOP;

    setInFib(false);
    setInTum(false);
    setInNormVes(false);
    setInTumVes(false);
    setInHypNec(false);
    setInMitCat(false);
    setInApop(false);

    ST_TIMER = IN_TIMER;

    if(ST_TUM){
        if(PAR_TUM_GROWTH){
            ST_G1 = ST_TIMER      <  PAR_LIM_G1S;
            ST_S  = PAR_LIM_G1S   <= ST_TIMER && ST_TIMER < PAR_LIM_SG2;
            ST_G2 = PAR_LIM_SG2   <= ST_TIMER && ST_TIMER < PAR_LIM_G2M;
            ST_M  = PAR_LIM_G2M   <= ST_TIMER && ST_TIMER < PAR_DOUB_TIME;
            ST_G0 = PAR_DOUB_TIME <= ST_TIMER;
        }
        else{
            ST_G1 = ST_TIMER == -1.0;
            ST_S  = ST_TIMER == -2.0;
            ST_G2 = ST_TIMER == -3.0;
            ST_M  = ST_TIMER == -4.0;
            ST_G0 = ST_TIMER == -5.0;
        }
    }

    else{
        ST_G1 = false;
        ST_S  = false;
        ST_G2 = false;
        ST_M  = false;
        ST_G0 = false;
    }

    if(ST_FIB){
        ST_ALPHA = PAR_ALPHA_FIB;
        ST_BETA  = PAR_BETA_FIB;
    }
    if(ST_TUM){
        if(ST_G1){
            ST_ALPHA = PAR_ALPHA_G1;
            ST_BETA  = PAR_BETA_G1;
        }
        if(ST_S){
            ST_ALPHA = PAR_ALPHA_S;
            ST_BETA  = PAR_BETA_S;
        }
        if(ST_G2){
            ST_ALPHA = PAR_ALPHA_G2;
            ST_BETA  = PAR_BETA_G2;
        }
        if(ST_M){
            ST_ALPHA = PAR_ALPHA_M;
            ST_BETA  = PAR_BETA_M;
        }
        if(ST_G0){
            ST_ALPHA = PAR_ALPHA_G0;
            ST_BETA  = PAR_BETA_G0;
        }
    }
    if(ST_NORM_VES){
        ST_ALPHA = PAR_ALPHA_NORM_VES;
        ST_BETA  = PAR_BETA_NORM_VES;
    }
    if(ST_TUM_VES){
        ST_ALPHA = PAR_ALPHA_TUM_VES;
        ST_BETA  = PAR_BETA_TUM_VES;
    }
    if(ST_DEAD){
        ST_ALPHA = 0.0;
        ST_BETA  = 0.0;
    }

    ST_ACC_DOSE = 0.0;

    ST_PO2  = IN_PO2;
    ST_VEGF = IN_VEGF;

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model terminalModel method.
------------------------------------------------------------------------------*/

int Cell::terminateModel(){
    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model updateModel method.
 *
 * Inputs:
 *  - currentTime: simulation current time (h),
 *  - DT: simulation timestep (h).
------------------------------------------------------------------------------*/

int Cell::updateModel(const double currentTime, const double DT){
    ST_PO2  = IN_PO2;
    ST_VEGF = IN_VEGF;

    if(PAR_OXY && !ST_DEAD){
        calcHypNec();
    }

    if(PAR_RES && ST_FIB){
        calcFibProlif(DT);
    }

    if(PAR_TUM_GROWTH && ST_TUM){
        calcTumGrowth(DT);
    }

    if(PAR_ANG){
        if(ST_NORM_VES){
            calcNormVesProlif(DT);
        }
        else if(ST_TUM_VES){
            calcTumVesProlif(DT);
        }
    }

    if(m_treatment){
        if(currentTime <= m_treatment->getDuration() &&
                fmod(currentTime, m_treatment->getInterval()) == DT){
            int i(currentTime / m_treatment->getInterval());
            if((m_treatment->getSchedule()).at(i)){
                ST_ACC_DOSE += m_treatment->getFraction();
                if(!ST_DEAD){
                    calcRespToIrr();
                }
            }
        }
    }

    if(IN_HYP_NEC){
        ST_FIB      = false;
        ST_TUM      = false;
        ST_NORM_VES = false;
        ST_TUM_VES  = false;
        ST_VES      = false;
        ST_HYP_NEC  = true;
        ST_MIT_CAT  = false;
        ST_APOP     = false;
        ST_DEAD     = true;

        ST_DAM    = false;
        ST_TIMER  = 0.0;
        ST_ARREST = 0.0;
    }

    else if(IN_MIT_CAT){
        ST_FIB      = false;
        ST_TUM      = false;
        ST_NORM_VES = false;
        ST_TUM_VES  = false;
        ST_VES      = false;
        ST_HYP_NEC  = false;
        ST_MIT_CAT  = true;
        ST_APOP     = false;
        ST_DEAD     = true;

        ST_DAM    = false;
        ST_TIMER  = 0.0;
        ST_ARREST = 0.0;
    }

    else if(IN_APOP){
        ST_FIB      = false;
        ST_TUM      = false;
        ST_NORM_VES = false;
        ST_TUM_VES  = false;
        ST_VES      = false;
        ST_HYP_NEC  = false;
        ST_MIT_CAT  = false;
        ST_APOP     = true;
        ST_DEAD     = true;

        ST_DAM    = false;
        ST_TIMER  = 0.0;
        ST_ARREST = 0.0;
    }

    else if(IN_NORM_VES){
        ST_FIB      = false;
        ST_TUM      = false;
        ST_NORM_VES = true;
        ST_TUM_VES  = false;
        ST_VES      = true;
        ST_HYP_NEC  = false;
        ST_MIT_CAT  = false;
        ST_APOP     = false;
        ST_DEAD     = false;

        ST_DAM    = false;
        ST_TIMER  = 0.0;
        ST_ARREST = 0.0;
    }

    else if(IN_TUM_VES){
        ST_FIB      = false;
        ST_TUM      = false;
        ST_NORM_VES = false;
        ST_TUM_VES  = true;
        ST_VES      = true;
        ST_HYP_NEC  = false;
        ST_MIT_CAT  = false;
        ST_APOP     = false;
        ST_DEAD     = false;

        ST_DAM    = false;
        ST_TIMER  = 0.0;
        ST_ARREST = 0.0;
    }

    else if(IN_FIB){
        ST_FIB      = true;
        ST_TUM      = false;
        ST_NORM_VES = false;
        ST_TUM_VES  = false;
        ST_VES      = false;
        ST_HYP_NEC  = false;
        ST_MIT_CAT  = false;
        ST_APOP     = false;
        ST_DEAD     = false;

        ST_DAM    = false;
        ST_TIMER  = 0.0;
        ST_ARREST = 0.0;
    }

    else if(IN_TUM){
        ST_FIB      = false;
        ST_TUM      = true;
        ST_NORM_VES = false;
        ST_TUM_VES  = false;
        ST_VES      = false;
        ST_HYP_NEC  = false;
        ST_MIT_CAT  = false;
        ST_APOP     = false;
        ST_DEAD     = false;

        ST_DAM    = false;
        ST_TIMER  = 0.0;
        ST_ARREST = 0.0;
        ST_G1     = true;
    }

    IN_FIB      = false;
    IN_TUM      = false;
    IN_NORM_VES = false;
    IN_TUM_VES  = false;
    IN_HYP_NEC  = false;
    IN_MIT_CAT  = false;
    IN_APOP     = false;

    if(!ST_TUM){
        ST_G1 = false;
        ST_S  = false;
        ST_G2 = false;
        ST_M  = false;
        ST_G0 = false;
    }

    if(ST_FIB){
        ST_ALPHA = PAR_ALPHA_FIB;
        ST_BETA  = PAR_BETA_FIB;
    }
    else if(ST_TUM){
        if(ST_G1){
            ST_ALPHA = PAR_ALPHA_G1;
            ST_BETA  = PAR_BETA_G1;
        }
        else if(ST_S){
            ST_ALPHA = PAR_ALPHA_S;
            ST_BETA  = PAR_BETA_S;
        }
        else if(ST_G2){
            ST_ALPHA = PAR_ALPHA_G2;
            ST_BETA  = PAR_BETA_G2;
        }
        else if(ST_M){
            ST_ALPHA = PAR_ALPHA_M;
            ST_BETA  = PAR_BETA_M;
        }
        else if(ST_G0){
            ST_ALPHA = PAR_ALPHA_G0;
            ST_BETA  = PAR_BETA_G0;
        }
    }
    else if(ST_NORM_VES){
        ST_ALPHA = PAR_ALPHA_NORM_VES;
        ST_BETA  = PAR_BETA_NORM_VES;
    }
    else if(ST_TUM_VES){
        ST_ALPHA = PAR_ALPHA_TUM_VES;
        ST_BETA  = PAR_BETA_TUM_VES;
    }
    else if(ST_DEAD){
        ST_ALPHA = 0.0;
        ST_BETA  = 0.0;
    }

    return 0;
}


/*------------------------------------------------------------------------------
 * This function adds a cell to the edge of the current one.
 *
 * Inputs:
 *  - cell: pointer to the cell to be added to the edge of the current one.
------------------------------------------------------------------------------*/

void Cell::addToEdge(Cell *const cell){
    m_edge->push_back(cell);
}


/*------------------------------------------------------------------------------
 * This function handles the progress in the cycle and potential division of
 * a healthy cell.
 *
 * Inputs:
 *  - DT: simulation timestep (h).
------------------------------------------------------------------------------*/

void Cell::calcFibProlif(double DT){
    if(ST_TIMER >= PAR_FIB_DOUB_TIME){
        Cell *newFib(0);
        newFib = searchSpaceForFib();
        if(newFib){
            ST_TIMER = 0.0;
            if(ST_DAM){
                IN_MIT_CAT = true;
            }
            else{
                ST_TIMER = 0.0;
                newFib->IN_FIB = true;
            }
        }
    }

    else{
        ST_TIMER += DT;
    }
}


/*------------------------------------------------------------------------------
 * This function handles the potential hypoxic necrosis of a cell.
 *
 * Inputs:
 *  - DT: simulation timestep (h).
------------------------------------------------------------------------------*/

void Cell::calcHypNec(){
    if(ST_PO2 < PAR_HYP_NEC_THRES){
        IN_HYP_NEC = true;
    }
}


/*------------------------------------------------------------------------------
 * This function handles the progress in the cycle and potential division of
 * a pre-existing endothelial cell.
 *
 * Inputs:
 *  - DT: simulation timestep (h).
------------------------------------------------------------------------------*/

void Cell::calcNormVesProlif(double DT){
    if(ST_TIMER >= PAR_ANG_TIME){
        if(ST_VEGF >= PAR_VEGF_THRES){
            Cell *newVes(0);
            newVes = searchSpaceForTumVes();
            if(newVes){
                if(ST_DAM){
                    IN_MIT_CAT = true;
                }
                else{
                    ST_TIMER = 0.0;
                    newVes->IN_TUM_VES = true;
                }
            }
        }
    }

    else{
        ST_TIMER += DT;
    }
}


/*------------------------------------------------------------------------------
 * This function calculates the Oxygen Enhancement Ratio of the response to
 * irradiation of a cell.
------------------------------------------------------------------------------*/

double Cell::calcOER() const{
    if(PAR_OXY){
        return (PAR_M * ST_PO2 + PAR_K) / (ST_PO2 + PAR_K); //mmHg
    }
    else{
        return PAR_M;
    }
}


/*------------------------------------------------------------------------------
 * This function handles the response to irradiation of a cell.
------------------------------------------------------------------------------*/

void Cell::calcRespToIrr(){
    ST_ARREST = PAR_ARREST_TIME;

    double p(double(rand()) / double(RAND_MAX));
    if(calcSF() < p){
        if(m_treatment->getFraction() < PAR_DOSE_THRES){
            if(PAR_TUM_GROWTH){
                ST_DAM = true;
            }
            else{
                IN_APOP = true;
            }
        }
        else{
            IN_APOP = true;
        }
    }
}


/*------------------------------------------------------------------------------
 * This function calculates the irradiation survival fraction of a cell.
------------------------------------------------------------------------------*/

double Cell::calcSF() const{
    double fraction, OER, SF;

    fraction = m_treatment->getFraction();
    OER = calcOER();
    SF = exp(-ST_ALPHA / PAR_M * fraction * OER -
             ST_BETA / (PAR_M * PAR_M) * fraction * fraction * OER * OER);
    return SF;
}


/*------------------------------------------------------------------------------
 * This function handles the progress in the cycle and potential mitotic death
 * or division of a tumour cell.
 *
 * Inputs:
 *  - DT: simulation timestep (h).
------------------------------------------------------------------------------*/

void Cell::calcTumGrowth(double DT){
    if(ST_G1 && ST_TIMER >= PAR_LIM_G1S){
        if(ST_ARREST > 0.0){
            ST_ARREST -= DT;
        }
        else{
            ST_TIMER += DT;
            ST_G1 = false;
            ST_S  = true;
        }
    }

    else if(ST_S && ST_TIMER >= PAR_LIM_SG2){
        ST_TIMER += DT;
        ST_S  = false;
        ST_G2 = true;
    }

    else if(ST_G2 && ST_TIMER >= PAR_LIM_G2M){
        if(ST_ARREST > 0.0){
            ST_ARREST -= DT;
        }
        else{
            ST_TIMER += DT;
            ST_G2 = false;
            ST_M  = true;
        }
    }

    else if(ST_TIMER >= PAR_DOUB_TIME){
        Cell *newTumCell(0);
        newTumCell = searchSpaceForTum();
        if(newTumCell){
            if(ST_DAM){
                IN_MIT_CAT = true;
            }
            else{
                ST_TIMER = 0.0;
                ST_G1 = true;
                ST_M  = false;
                ST_G0 = false;
                newTumCell->IN_TUM = true;
            }
        }
        else{
            ST_M  = false;
            ST_G0 = true;
        }
    }

    else{
        ST_TIMER += DT;
    }
}


/*------------------------------------------------------------------------------
 * This function handles the progress in the cycle and potential division of
 * a neo-created endothelial cell.
 *
 * Inputs:
 *  - DT: simulation timestep (h).
------------------------------------------------------------------------------*/

void Cell::calcTumVesProlif(double DT){
    if(ST_TIMER >= PAR_ANG_TIME){
        if(ST_VEGF >= PAR_VEGF_THRES){
            Cell *newVes(0);
            newVes = searchSpaceForTumVes();
            if(newVes){
                if(ST_DAM){
                    IN_MIT_CAT = true;
                }
                else{
                    ST_TIMER = 0.0;
                    newVes->IN_TUM_VES = true;
                }
            }
        }
    }
    else{
        ST_TIMER += DT;
    }
}


/*------------------------------------------------------------------------------
 * This function gets the healthy cell state.
------------------------------------------------------------------------------*/

bool Cell::getFib() const{
    return ST_FIB;
}


/*------------------------------------------------------------------------------
 * This function gets accumulated dose.
------------------------------------------------------------------------------*/

double Cell::getAccDose() const{
    return ST_ACC_DOSE;
}


/*------------------------------------------------------------------------------
 * This function gets the healthy cell state.
------------------------------------------------------------------------------*/

bool Cell::getDead() const{
    return ST_DEAD;
}


/*------------------------------------------------------------------------------
 * This function gets the edge of the current cell.
------------------------------------------------------------------------------*/

vector<Cell *> *Cell::getEdge() const{
    return m_edge;
}


/*------------------------------------------------------------------------------
 * This function gets the G0 phase state.
------------------------------------------------------------------------------*/

bool Cell::getG0() const{
    return ST_G0;
}


/*------------------------------------------------------------------------------
 * This function gets the G1 phase state.
------------------------------------------------------------------------------*/

bool Cell::getG1() const{
    return ST_G1;
}


/*------------------------------------------------------------------------------
 * This function gets the G2 phase state.
------------------------------------------------------------------------------*/

bool Cell::getG2() const{
    return ST_G2;
}


/*------------------------------------------------------------------------------
 * This function gets the hypoxic necrotic cell state.
------------------------------------------------------------------------------*/

bool Cell::getHypNec() const{
    return ST_HYP_NEC;
}


/*------------------------------------------------------------------------------
 * This function gets the M phase state.
------------------------------------------------------------------------------*/

bool Cell::getM() const{
    return ST_M;
}


/*------------------------------------------------------------------------------
 * This function gets the mitotic dead cell state.
------------------------------------------------------------------------------*/

bool Cell::getMitCat() const{
    return ST_MIT_CAT;
}


/*------------------------------------------------------------------------------
 * This function gets the pre-existing endothelial cell state.
------------------------------------------------------------------------------*/

bool Cell::getNormVes() const{
    return ST_NORM_VES;
}


/*------------------------------------------------------------------------------
 * This function gets state output.
------------------------------------------------------------------------------*/

int Cell::getOutState() const{
    return OUT_STATE;
}


/*------------------------------------------------------------------------------
 * This function gets the S phase state.
------------------------------------------------------------------------------*/

bool Cell::getS() const{
    return ST_S;
}


/*------------------------------------------------------------------------------
 * This function gets the tumour cell state.
------------------------------------------------------------------------------*/

bool Cell::getTum() const{
    return ST_TUM;
}


/*------------------------------------------------------------------------------
 * This function gets the neo-created endothelial cell state.
------------------------------------------------------------------------------*/

bool Cell::getTumVes() const{
    return ST_TUM_VES;
}


/*------------------------------------------------------------------------------
 * This function gets the damaged tumour cell state.
------------------------------------------------------------------------------*/

bool Cell::getTumDam() const{
    return ST_TUM && ST_DAM;
}


/*------------------------------------------------------------------------------
 * This function gets the not damaged tumour cell state.
------------------------------------------------------------------------------*/

bool Cell::getTumNotDam() const{
    return ST_TUM && !ST_DAM;
}


/*------------------------------------------------------------------------------
 * This function gets the endothelial cell state.
------------------------------------------------------------------------------*/

bool Cell::getVes() const{
    return ST_VES;
}


/*------------------------------------------------------------------------------
 * This function searches an available space for an initial tumour cell.
------------------------------------------------------------------------------*/

Cell *Cell::searchInitSpaceForTum() const{
    int edgeSize, m;

    edgeSize = m_edge->size();
    m = rand() % edgeSize;
    for(int n(0); n < edgeSize; n++){
        if(!m_edge->at(m)->IN_TUM && !m_edge->at(m)->IN_NORM_VES &&
                !m_edge->at(m)->IN_TUM_VES){
            return m_edge->at(m);
        }
        m++;
        if(m == edgeSize){
            m = 0;
        }
    }
    return 0;
}


/*------------------------------------------------------------------------------
 * This function searches an available space for a healthy cell.
------------------------------------------------------------------------------*/

Cell *Cell::searchSpaceForFib() const{
    int edgeSize, m;

    edgeSize = m_edge->size();
    m = rand() % edgeSize;
    for(int n(0); n < edgeSize; n++){
        if((m_edge->at(m)->getOutState() == 6 ||
            m_edge->at(m)->getOutState() == 7 ||
            m_edge->at(m)->getOutState() == 8) && !m_edge->at(m)->IN_FIB &&
                !m_edge->at(m)->IN_TUM && !m_edge->at(m)->IN_NORM_VES &&
                !m_edge->at(m)->IN_TUM_VES){
            return m_edge->at(m);
        }
        m++;
        if(m == edgeSize){
            m = 0;
        }
    }
    return 0;
}


/*------------------------------------------------------------------------------
 * This function searches an available space for a tumour cell.
------------------------------------------------------------------------------*/

Cell *Cell::searchSpaceForTum() const{
    int edgeSize, m;

    edgeSize = m_edge->size();
    m = rand() % edgeSize;
    for(int n(0); n < edgeSize; n++){
        if((m_edge->at(m)->getOutState() == 1 ||
            m_edge->at(m)->getOutState() == 7 ||
            m_edge->at(m)->getOutState() == 8) && !m_edge->at(m)->IN_FIB &&
                !m_edge->at(m)->IN_TUM && !m_edge->at(m)->IN_NORM_VES &&
                !m_edge->at(m)->IN_TUM_VES){
            return m_edge->at(m);
        }
        m++;
        if(m == edgeSize){
            m = 0;
        }
    }
    return 0;
}


/*------------------------------------------------------------------------------
 * This function searches an available space for a neo-created endothelial cell.
------------------------------------------------------------------------------*/

Cell *Cell::searchSpaceForTumVes() const{
    int edgeSize;
    double m(0.0);
    Cell *newTumVes(0), *tempCell(0);

    edgeSize = m_edge->size();

    for(int n(0); n < edgeSize; n++){
        tempCell = m_edge->at(n);
        if(tempCell->getOutState() != 4 && tempCell->getOutState() != 5 &&
                tempCell->getOutState() != 6 && !tempCell->IN_FIB &&
                !tempCell->IN_TUM && !tempCell->IN_NORM_VES &&
                !tempCell->IN_TUM_VES){
            if(tempCell->ST_VEGF > m){
                m = tempCell->ST_VEGF;
                newTumVes = tempCell;
            }
        }
    }
    return newTumVes;
}


/*------------------------------------------------------------------------------
 * This function sets the healthy cell input.
 *
 * Inputs:
 *  - input: healthy cell input.
------------------------------------------------------------------------------*/

void Cell::setInFib(const bool input){
    IN_FIB = input;
}


/*------------------------------------------------------------------------------
 * This function sets the apoptotic cell input.
 *
 * Inputs:
 *  - input: apoptotic cell input.
------------------------------------------------------------------------------*/

void Cell::setInApop(const bool input){
    IN_APOP = input;
}


/*------------------------------------------------------------------------------
 * This function sets the hypoxic necrotic cell input.
 *
 * Inputs:
 *  - input: hypoxic necrotic cell input.
------------------------------------------------------------------------------*/

void Cell::setInHypNec(const bool input){
    IN_HYP_NEC = input;
}


/*------------------------------------------------------------------------------
 * This function sets the mitotic dead cell input.
 *
 * Inputs:
 *  - input: mitotic dead cell input.
------------------------------------------------------------------------------*/

void Cell::setInMitCat(const bool input){
    IN_MIT_CAT = input;
}


/*------------------------------------------------------------------------------
 * This function sets the pre-existing endothelial cell input.
 *
 * Inputs:
 *  - input: pre-existing endothelial cell input.
------------------------------------------------------------------------------*/

void Cell::setInNormVes(const bool input){
    IN_NORM_VES = input;
}


/*------------------------------------------------------------------------------
 * This function sets the pO2 input.
 *
 * Inputs:
 *  - input: pO2 input (mmHg).
------------------------------------------------------------------------------*/

void Cell::setInPO2(const double input){
    IN_PO2 = input;
}


/*------------------------------------------------------------------------------
 * This function sets the timer input.
 *
 * Inputs:
 *  - input: timer input (h).
------------------------------------------------------------------------------*/

void Cell::setInTimer(const double input){
    IN_TIMER = input;
}


/*------------------------------------------------------------------------------
 * This function sets the tumour cell input.
 *
 * Inputs:
 *  - input: tumour cell input.
------------------------------------------------------------------------------*/

void Cell::setInTum(const bool input){
    IN_TUM = input;
}


/*------------------------------------------------------------------------------
 * This function sets the neo-created endothelial cell input.
 *
 * Inputs:
 *  - input: neo-created endothelial cell input.
------------------------------------------------------------------------------*/

void Cell::setInTumVes(const bool input){
    IN_TUM_VES = input;
}


/*------------------------------------------------------------------------------
 * This function sets the VEGF concentration input.
 *
 * Inputs:
 *  - input: VEGF concentration input (mol/um^3).
------------------------------------------------------------------------------*/

void Cell::setInVegf(const double input){
    IN_VEGF = input;
}
