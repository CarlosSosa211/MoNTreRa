/**
 * @file tissue.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.19.17
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "tissue.hpp"

using namespace std;

/*------------------------------------------------------------------------------
 * Constructor of the class Tissue.
 *
 * Inputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue (um),
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - tumGrowth: tumour growth,
 *  - doubTime: duration of the cycle of tumour cells (h),
 *  - edgeOrder: order of the edge of cells,
 *  - cycDur: vector containing the duration fractions of every phase of the
 *  cell cycle,
 *  - cycDistrib: vector containing the initial distribution of tumour cells in
 *  the cycle,
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
 *  - treatment: pointer to the treatment,
 *  - doseThres: dose threshold to provoke instantaneous death by apoptosis
 *  (Gy),
 *  - arrestTime: radiation-induced arrest time (h),
 *  - oxy: integer indicating the oxygenation scenario (1, space-and-time
 *  dependent),
 *  - hypNecThres: pO2 hypoxic necrosis thresthold (mmHg).
------------------------------------------------------------------------------*/

Tissue::Tissue(const int nrow, const int ncol, const int nlayer,
               const double cellSize, const vector<bool> &inTum,
               const vector<bool> &inVes, const bool tumGrowth,
               const double doubTime, const int edgeOrder,
               vector<double> cycDur, vector<double> cycDistrib,
               const bool res, const double fibDoubTime, const bool ang,
               const double angTime, const double vegfThres,
               vector<double> alpha, vector<double> beta,
               Treatment *const treatment, const double doseThres,
               const double arrestTime, const int oxy,
               const double hypNecThres) :
    Model(TISSUE_NUM_IN_B, TISSUE_NUM_IN_I, TISSUE_NUM_IN_D, TISSUE_NUM_ST_B,
          TISSUE_NUM_ST_I, TISSUE_NUM_ST_D, TISSUE_NUM_OUT_B, TISSUE_NUM_OUT_I,
          TISSUE_NUM_OUT_D, TISSUE_NUM_PAR_B, TISSUE_NUM_PAR_I,
          TISSUE_NUM_PAR_D, nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;
    m_cellSize = cellSize * 1e-3; //(mm)

    m_treatment = treatment;

    int k(0);
    double inputTimer, n;

    srand(time(NULL));

    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j++){
                m_comp->at(k) = new Cell(i, j, l, tumGrowth, doubTime, cycDur,
                                         res, fibDoubTime, ang, angTime,
                                         vegfThres, alpha, beta, doseThres,
                                         arrestTime, oxy, hypNecThres, this);

                ((Cell *)m_comp->at(k))->setInTum(inTum.at(k));

                if(inTum.at(k) && inVes.at(k)){
                    cout << "Conflict between initial data. Cell "<< k <<
                            " is both tumor and vessel" << endl;
                    break;
                }

                else{
                    ((Cell *)m_comp->at(k))->setInNormVes(inVes.at(k));
                }
                k++;
            }
        }
    }

    int incol, iincol, lnrowNcol, llnrowNcol;
    for(int l(0); l < m_nlayer; l++){
        lnrowNcol = l * m_nrow * m_ncol;
        for(int i(0); i < m_nrow; i++){
            incol = i * m_ncol;
            for(int j(0); j < m_ncol; j++){
                for(int ll(-edgeOrder); ll <= edgeOrder; ll++){
                    llnrowNcol = lnrowNcol + ll * m_nrow * m_ncol;
                    for(int ii(-edgeOrder); ii <= edgeOrder; ii++){
                        iincol = llnrowNcol + incol + ii * m_ncol;
                        for(int jj(-edgeOrder); jj <= edgeOrder; jj++){
                            if(ii != 0 || jj != 0 || ll != 0){
                                if(l + ll >= 0 && l + ll < m_nlayer &&
                                        i + ii >= 0 && i + ii < m_nrow &&
                                        j + jj >= 0 && j + jj < m_ncol){
                                    ((Cell *)m_comp->at(lnrowNcol + incol + j))
                                            ->addToEdge((Cell *)m_comp->
                                                        at(iincol + j + jj));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for(int k(0); k < m_numComp; k++){
        if(inTum.at(k)){
            if(tumGrowth){
                if(((Cell *)m_comp->at(k))->searchInitSpaceForTum()){
                    n = double(rand()) / double(RAND_MAX);
                    if(n < cycDistrib.at(0)){
                        inputTimer = rand() % int(cycDur.at(0) * doubTime);
                    }
                    else if(n < cycDistrib.at(0) + cycDistrib.at(1)){
                        inputTimer = cycDur.at(0) * doubTime +
                                rand() % int(cycDur.at(1) * doubTime);
                    }
                    else if(n < cycDistrib.at(0) + cycDistrib.at(1) +
                            cycDistrib.at(2)){
                        inputTimer = (cycDur.at(0) + cycDur.at(1)) * doubTime +
                                rand() % int(cycDur.at(2) * doubTime);
                    }
                    else{
                        inputTimer = (cycDur.at(0) + cycDur.at(1) +
                                      cycDur.at(2)) * doubTime + rand() %
                                int(cycDur.at(3) * doubTime);
                    }
                }
                else{
                    inputTimer = doubTime;
                }
            }

            else{
                if(((Cell *)m_comp->at(k))->searchInitSpaceForTum()){
                    n = double(rand()) / double(RAND_MAX);
                    if(n < cycDistrib.at(0)){
                        inputTimer = -1.0;
                    }
                    else if(n < cycDistrib.at(0) + cycDistrib.at(1)){
                        inputTimer = -2.0;
                    }
                    else if(n < cycDistrib.at(0) + cycDistrib.at(1) +
                            cycDistrib.at(2)){
                        inputTimer = -3.0;
                    }
                    else{
                        inputTimer = -4.0;
                    }
                }
                else{
                    inputTimer = -5.0;
                }
            }
        }

        else if(inVes.at(k)){
            if(ang){
                n = double(rand()) / double(RAND_MAX);
                inputTimer = n * angTime;
            }
            else{
                inputTimer = 0.0;
            }
        }

        else{
            if(1){
                n = double(rand()) / double(RAND_MAX);
                inputTimer = n * fibDoubTime;
            }
            else{
                inputTimer = 0.0;
            }
        }

        ((Cell *)m_comp->at(k))->setInTimer(inputTimer);
    }
}

/*------------------------------------------------------------------------------
 * Destructor of the class Tissue.
------------------------------------------------------------------------------*/

Tissue::~Tissue(){
    for(int k(0); k < m_numComp; k++){
        delete m_comp->at(k);
    }
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model calcModelOut method.
------------------------------------------------------------------------------*/

int Tissue::calcModelOut(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->calcModelOut();
    }

    OUT_TUM_DENS           = ST_TUM_DENS;
    OUT_END_TREAT_TUM_DENS = ST_END_TREAT_TUM_DENS;
    OUT_3MON_TUM_DENS      = ST_3MON_TUM_DENS;
    OUT_INT_TUM_DENS       = ST_INT_TUM_DENS;

    OUT_REC          = ST_REC;
    OUT_REC_TUM_DENS = ST_REC_TUM_DENS;
    OUT_REC_TIME     = ST_REC_TIME;

    if(PAR_INIT_TUM_DENS){
        OUT_KILLED_CELLS = (PAR_INIT_TUM_DENS - ST_TUM_DENS) /
                PAR_INIT_TUM_DENS * 100.0;
    }

    OUT_VES_DENS      = ST_VES_DENS;
    OUT_NORM_VES_DENS = ST_NORM_VES_DENS;
    OUT_TUM_VES_DENS  = ST_TUM_VES_DENS;

    OUT_DEAD_DENS = ST_DEAD_DENS;

    const int numTum(getNumTum());

    if(numTum){
        const double numTum100(100.0 / double(numTum));
        OUT_G1_DENS = numTum100 * double(getNumG1());
        OUT_S_DENS  = numTum100 * double(getNumS());
        OUT_G2_DENS = numTum100 * double(getNumG2());
        OUT_M_DENS  = numTum100 * double(getNumM());
        OUT_G0_DENS = numTum100 * double(getNumG0());
    }
    else{
        OUT_G1_DENS = 0.0;
        OUT_S_DENS  = 0.0;
        OUT_G2_DENS = 0.0;
        OUT_M_DENS  = 0.0;
        OUT_G0_DENS = 0.0;
    }

    OUT_50_KILLED  = ST_50_KILLED;
    OUT_80_KILLED  = ST_80_KILLED;;
    OUT_90_KILLED  = ST_90_KILLED;;
    OUT_95_KILLED  = ST_95_KILLED;;
    OUT_99_KILLED  = ST_99_KILLED;;
    OUT_999_KILLED = ST_999_KILLED;;

    OUT_TIME_TO_50  = ST_TIME_TO_50;
    OUT_TIME_TO_80  = ST_TIME_TO_80;
    OUT_TIME_TO_90  = ST_TIME_TO_90;
    OUT_TIME_TO_95  = ST_TIME_TO_95;
    OUT_TIME_TO_99  = ST_TIME_TO_99;
    OUT_TIME_TO_999 = ST_TIME_TO_999;

    OUT_DOSE_TO_50  = ST_DOSE_TO_50;
    OUT_DOSE_TO_80  = ST_DOSE_TO_80;
    OUT_DOSE_TO_90  = ST_DOSE_TO_90;
    OUT_DOSE_TO_95  = ST_DOSE_TO_95;
    OUT_DOSE_TO_99  = ST_DOSE_TO_99;
    OUT_DOSE_TO_999 = ST_DOSE_TO_999;

    OUT_CONTROLLED      = ST_CONTROLLED;
    OUT_DOSE_TO_CONTROL = ST_DOSE_TO_CONTROL;

    if(m_nlayer == 1){
        OUT_TUM_VOL = 4.0 / (3.0 * sqrt(M_PI)) * pow(numTum * m_cellSize *
                                                     m_cellSize, 1.5);
    }
    else{
        OUT_TUM_VOL = numTum * m_cellSize * m_cellSize * m_cellSize;
    }

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model initModel method.
------------------------------------------------------------------------------*/

int Tissue::initModel(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->initModel();
    }

    double _numComp100(1.0 / double(m_numComp) * 100.0);

    PAR_INIT_TUM_DENS = double(getNumTum()) * _numComp100;
    PAR_INIT_VES_DENS = double(getNumVes()) * _numComp100;

    ST_TUM_DENS           = PAR_INIT_TUM_DENS;
    ST_PREV_TUM_DENS      = PAR_INIT_TUM_DENS;
    ST_INT_TUM_DENS       = 0.0;
    ST_END_TREAT_TUM_DENS = 0.0;
    ST_3MON_TUM_DENS      = 0.0;

    ST_VES_DENS      = double(getNumVes()) * _numComp100;
    ST_NORM_VES_DENS = double(getNumNormVes()) * _numComp100;
    ST_TUM_VES_DENS  = double(getNumTumVes()) * _numComp100;

    ST_DEAD_DENS  = double(getNumDead()) * _numComp100;

    ST_REC          = false;
    ST_COUNT_REC    = 0;
    ST_REC_TUM_DENS = 0.0;
    ST_REC_TIME     = 0.0;

    ST_50_KILLED  = false;
    ST_80_KILLED  = false;
    ST_90_KILLED  = false;
    ST_95_KILLED  = false;
    ST_99_KILLED  = false;
    ST_999_KILLED = false;

    ST_TIME_TO_50  = 0.0;
    ST_TIME_TO_80  = 0.0;
    ST_TIME_TO_90  = 0.0;
    ST_TIME_TO_95  = 0.0;
    ST_TIME_TO_99  = 0.0;
    ST_TIME_TO_999 = 0.0;

    ST_DOSE_TO_50  = 0.0;
    ST_DOSE_TO_80  = 0.0;
    ST_DOSE_TO_90  = 0.0;
    ST_DOSE_TO_95  = 0.0;
    ST_DOSE_TO_99  = 0.0;
    ST_DOSE_TO_999 = 0.0;

    ST_CONTROLLED = false;

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model terminalModel method.
------------------------------------------------------------------------------*/

int Tissue::terminateModel(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->terminateModel();
    }
    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model updateModel method.
 *
 * Inputs:
 *  - currentTime: simulation current time (h),
 *  - DT: simulation timestep (h).
------------------------------------------------------------------------------*/

int Tissue::updateModel(const double currentTime, const double DT){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->updateModel(currentTime, DT);
    }

    double _numComp100(1.0 / double(m_numComp) * 100.0);
    ST_PREV_TUM_DENS = ST_TUM_DENS;
    ST_TUM_DENS      = double(getNumTum()) * _numComp100;
    ST_INT_TUM_DENS += 0.5 * DT * (ST_PREV_TUM_DENS + ST_TUM_DENS);

    ST_VES_DENS      = double(getNumVes()) * _numComp100;
    ST_NORM_VES_DENS = double(getNumNormVes()) * _numComp100;
    ST_TUM_VES_DENS  = double(getNumTumVes()) * _numComp100;

    ST_DEAD_DENS  = double(getNumDead()) * _numComp100;

    if(m_treatment){
        if(currentTime <= m_treatment->getDuration()){
            ST_END_TREAT_TUM_DENS = ST_TUM_DENS;
        }
        else if(!ST_REC && ST_TUM_DENS > ST_PREV_TUM_DENS){
            if(!ST_COUNT_REC){
                ST_REC_TUM_DENS = ST_TUM_DENS;
                ST_REC_TIME = currentTime;
            }
            ST_COUNT_REC ++;
            if(ST_COUNT_REC = 10){
                ST_REC = true;
            }
        }
        else{
            ST_COUNT_REC = 0;
        }
        if(currentTime <= 2160.0){
            ST_3MON_TUM_DENS = ST_TUM_DENS;
        }

        double tumSurv;
        tumSurv = ST_TUM_DENS / PAR_INIT_TUM_DENS;
        if(tumSurv < 0.5 && !ST_50_KILLED){
            ST_50_KILLED = true;
            ST_TIME_TO_50 = currentTime;
            ST_DOSE_TO_50 = ((Cell*)m_comp->at(0))->getAccDose();
        }
        if(tumSurv < 0.2 && !ST_80_KILLED){
            ST_80_KILLED = true;
            ST_TIME_TO_80 = currentTime;
            ST_DOSE_TO_80 = ((Cell*)m_comp->at(0))->getAccDose();
        }
        if(tumSurv < 0.1 && !ST_90_KILLED){
            ST_90_KILLED = true;
            ST_TIME_TO_90 = currentTime;
            ST_DOSE_TO_90 = ((Cell*)m_comp->at(0))->getAccDose();
        }
        if(tumSurv < 0.05 && !ST_95_KILLED){
            ST_95_KILLED = true;
            ST_TIME_TO_95 = currentTime;
            ST_DOSE_TO_95 = ((Cell*)m_comp->at(0))->getAccDose();
        }
        if(tumSurv < 0.01 && !ST_99_KILLED){
            ST_99_KILLED = true;
            ST_TIME_TO_99 = currentTime;
            ST_DOSE_TO_99 = ((Cell*)m_comp->at(0))->getAccDose();
        }
        if(tumSurv < 0.001 && !ST_999_KILLED){
            ST_999_KILLED = true;
            ST_TIME_TO_999 = currentTime;
            ST_DOSE_TO_999 = ((Cell*)m_comp->at(0))->getAccDose();
        }

        if(!getNumTumNotDam() && !ST_CONTROLLED){
            ST_CONTROLLED = true;
            ST_DOSE_TO_CONTROL = ((Cell*)m_comp->at(0))->getAccDose();
        }
    }
    return 0;
}


/*------------------------------------------------------------------------------
 * This function counts the healthy cells of the tissue.
 *
 * Outputs:
 *  - count: number of healthy cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumFib() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getFib()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the dead cells of the tissue.
 *
 * Outputs:
 *  - count: number of dead cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumDead() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getDead()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the G0 cells of the tissue.
 *
 * Outputs:
 *  - count: number of G0 cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumG0() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getG0()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the G1 cells of the tissue.
 *
 * Outputs:
 *  - count: number of G1 cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumG1() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getG1()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the G2 cells of the tissue.
 *
 * Outputs:
 *  - count: number of G2 cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumG2() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getG2()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the M cells of the tissue.
 *
 * Outputs:
 *  - count: number of M cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumM() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getM()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the pre-existing endothelial cells of the tissue.
 *
 * Outputs:
 *  - count: number of pre-existing endothelial cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumNormVes() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getNormVes()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the S cells of the tissue.
 *
 * Outputs:
 *  - count: number of S cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumS() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getS()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the tumour cells of the tissue.
 *
 * Outputs:
 *  - count: number of tumour cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumTum() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getTum()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the not damaged tumour cells of the tissue.
 *
 * Outputs:
 *  - count: number of not damaged tumour cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumTumNotDam() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getTumNotDam()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the neo-created endothelial cells of the tissue.
 *
 * Outputs:
 *  - count: number of neo-created endothelial cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumTumVes() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getTumVes()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function counts the endothelial cells of the tissue.
 *
 * Outputs:
 *  - count: number of endothelial cells of the tissue
------------------------------------------------------------------------------*/

int Tissue::getNumVes() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getVes()){
            count++;
        }
    }
    return count;
}


/*------------------------------------------------------------------------------
 * This function gets the treatment.
------------------------------------------------------------------------------*/

Treatment *Tissue::getTreatment() const{
    return m_treatment;
}


