/**
 * @file tissue.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "tissue.hpp"

using namespace std;

Tissue::Tissue(const int nrow, const int ncol,
               const int nlayer,
               Treatment *const treatment) :
    Model(0, 3, 23, 2, nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    m_treatment = treatment;

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new Cell(this);
        m_numOut += (m_comp->at(k))->getNumOut();
    }
}


Tissue::Tissue(const int nrow, const int ncol, const int nlayer,
               const string nFInTum, const string nFInVes,
               const double tumGrowth, const double doubTime,
               const int edgeOrder, vector<double> cycDur,
               vector<double> cycDistrib, const double res,
               const double fibDoubTime, const double ang,
               const double angTime, const double vegfThres,
               vector<double> alpha, vector<double> beta,
               const double doseThres, const double arrestTime, 
               Treatment *const treatment, const double hypNecThres) :
    Model(0, 3, 23, 2, nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    m_treatment = treatment;

    //Creation of the cells composing the tissue model
    int k(0);
    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j++){
                m_comp->at(k) = new Cell(i, j, l, tumGrowth, doubTime, cycDur,
                                         res, fibDoubTime, ang, angTime, vegfThres,
                                         alpha, beta, doseThres, arrestTime, 
                                         hypNecThres, this);
                m_numOut += (m_comp->at(k))->getNumOut();
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
                                if(l + ll >= 0 && l + ll < m_nlayer && i + ii >= 0
                                        && i + ii < m_nrow && j + jj >= 0 &&
                                        j + jj < m_ncol){
                                    ((Cell *)m_comp->at(lnrowNcol + incol + j))
                                            ->addToEdge((Cell *)m_comp->at(iincol + j + jj));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    double inputTimer, inputTum, inputVes, n;
    ifstream fInTum(nFInTum.c_str());
    ifstream fInVes(nFInVes.c_str());

    //Initialization of the cells state
    if(!fInTum.is_open()){
        /*cout << "An error occurred while opening initial tumor" <<
                "data file" << endl*/;
    }
    else if(!fInVes.is_open()){
        /*cout << "An error occurred while opening initial vessel" <<
                "data file" << endl*/;
    }
    else{
        srand(time(NULL));
        for(int k(0); k < m_numComp; k++){
            if(fInTum >> inputTum){
                ((Cell *)m_comp->at(k))->setInTum(inputTum);
            }
            else{
                /*cout << "Insufficient data in tumor file" << endl*/;
                break;
            }
            if(fInVes >> inputVes){
                if(inputTum && inputVes){
                    /*cout << "Conflict between initial data. Cell "<< k <<
                            " is both tumor and vessel" << endl*/;
                    break;
                }
                else{
                    ((Cell *)m_comp->at(k))->setInNormVes(inputVes);
                }
            }
            else{
                /*cout << "Insufficient data in vessel file" << endl*/;
                break;
            }

            if(inputTum){
                if(tumGrowth){
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
                        inputTimer = (cycDur.at(0) + cycDur.at(1) + cycDur.at(2)) *
                                doubTime + rand() % int(cycDur.at(3) * doubTime);
                    }
                }

                else{
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
            }

            else if(inputVes){
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
                    inputTimer = n * 500.0;
                }
                else{
                    inputTimer = 0.0;
                }
            }

            ((Cell *)m_comp->at(k))->setInTimer(inputTimer);
        }
        fInTum.close();
        fInVes.close();
    }
}


Tissue::Tissue(const int nrow, const int ncol, const int nlayer,
               const vector<bool> &inTum, const vector<bool> &inVes,
               const double tumGrowth, const double doubTime,
               const int edgeOrder, vector<double> cycDur,
               vector<double> cycDistrib, const double res,
               const double fibDoubTime, const double ang,
               const double angTime, const double vegfThres,
               vector<double> alpha, vector<double> beta,
               const double doseThres, const double arrestTime,
               Treatment *const treatment, const double hypNecThres) :
    Model(0, 3, 23, 2, nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    m_treatment = treatment;

    //Creation of the cells composing the tissue model
    int k(0);
    double inputTimer, n;

    srand(time(NULL));

    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j++){
                m_comp->at(k) = new Cell(i, j, l, tumGrowth, doubTime, cycDur,
                                         res, fibDoubTime, ang, angTime,
                                         vegfThres, alpha, beta, doseThres, 
                                         arrestTime, hypNecThres, this);
                m_numOut += (m_comp->at(k))->getNumOut();

                ((Cell *)m_comp->at(k))->setInTum(inTum.at(k));

                if(inTum.at(k) && inVes.at(k)){
                    /*cout << "Conflict between initial data. Cell "<< k <<
                            " is both tumor and vessel" << endl*/;
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
                                if(l + ll >= 0 && l + ll < m_nlayer && i + ii >= 0
                                        && i + ii < m_nrow && j + jj >= 0 &&
                                        j + jj < m_ncol){
                                    ((Cell *)m_comp->at(lnrowNcol + incol + j))
                                            ->addToEdge((Cell *)m_comp->at(iincol + j + jj));
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
                        inputTimer = (cycDur.at(0) + cycDur.at(1) + cycDur.at(2)) *
                                doubTime + rand() % int(cycDur.at(3) * doubTime);
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


Tissue::~Tissue(){
    for(int k(0); k < m_numComp; k++){
        delete m_comp->at(k);
    }
}


int Tissue::calcModelOut(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->calcModelOut();
    }

    OUT_TUM_DENS     = ST_TUM_DENS;
    OUT_INT_TUM_DENS = ST_INT_TUM_DENS;
    if(PAR_INIT_TUM_DENS){
        OUT_KILLED_CELLS = (PAR_INIT_TUM_DENS - ST_TUM_DENS) / PAR_INIT_TUM_DENS * 100.0;
    }
    OUT_VES_DENS      = double(getNumVes()) / double(m_numComp) * 100.0;
    OUT_NORM_VES_DENS = double(getNumNormVes()) / double(m_numComp) * 100.0;
    OUT_TUM_VES_DENS  = double(getNumTumVes()) / double(m_numComp) * 100.0;

    if(getNumTum()){
        double numTum100(100.0 / double(getNumTum()));
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

    OUT_TIME_TO_50  = m_timeNeeded[0];
    OUT_TIME_TO_80  = m_timeNeeded[1];
    OUT_TIME_TO_90  = m_timeNeeded[2];
    OUT_TIME_TO_95  = m_timeNeeded[3];
    OUT_TIME_TO_99  = m_timeNeeded[4];
    OUT_TIME_TO_999 = m_timeNeeded[5];

    OUT_DOSE_TO_50  = m_doseNeeded[0];
    OUT_DOSE_TO_80  = m_doseNeeded[1];
    OUT_DOSE_TO_90  = m_doseNeeded[2];
    OUT_DOSE_TO_95  = m_doseNeeded[3];
    OUT_DOSE_TO_99  = m_doseNeeded[4];
    OUT_DOSE_TO_999 = m_doseNeeded[5];

    return 0;
}


int Tissue::initModel(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->initModel();
    }

    PAR_INIT_TUM_DENS = double(getNumTum()) / double(m_numComp) * 100.0;
    PAR_INIT_VES_DENS = double(getNumVes()) / double(m_numComp) * 100.0;

    ST_TUM_DENS      = PAR_INIT_TUM_DENS;
    ST_PREV_TUM_DENS = PAR_INIT_TUM_DENS;
    ST_INT_TUM_DENS  = 0.0;

    for(int i(0); i < 6; i++){
        //m_doseNeeded[i] = m_treatment->getTotalDose() + m_treatment->getFraction();
        //m_timeNeeded[i] = m_treatment->getDuration() + m_treatment->getInterval();
        m_doseNeeded[i] = -1.0;
        m_timeNeeded[i] = -1.0;
    }

    /*cout << "Total number of cells = " << m_numComp << endl*/;
    /*cout << "Initial number of cells at G1 = " << getNumG1() << endl*/;
    /*cout << "Initial number of cells at S = " << getNumS() << endl*/;
    /*cout << "Initial number of cells at G2 = " << getNumG2() << endl*/;
    /*cout << "Initial number of cells at M = " << getNumM() << endl*/;
    /*cout << "Initial number of cells at G0 = " << getNumG0() << endl*/;
    /*cout << "Initial number of living cells = "
         << getNumFib() << endl*/;
    /*cout << "Initial number of tumor cells = " << getNumTum() << endl*/;
    /*cout << "Initial tumor density: " << PAR_INIT_TUM_DENS << "%"
         << endl*/;
    /*cout << "Initial number of vessels = " << getNumVes() << endl*/;
    /*cout << "Initial vascular density: " << PAR_INIT_VES_DENS << "%"
         << endl*/;
    /*cout << "Initial number of dead cells = " << getNumDead() << endl*/;
    /*cout << "---------------------------------------------" << endl*/;

    return 0;
}


int Tissue::terminateModel(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->terminateModel();
    }
    if(m_treatment){
        for(int i(0); i < 6; i++){
            if(m_doseNeeded[i] < 0){
                m_doseNeeded[i] = m_treatment->getTotalDose() + m_treatment->getFraction();
                m_timeNeeded[i] = m_treatment->getDuration() + m_treatment->getInterval();
            }
        }
    }
    /*cout << "---------------------------------------------" << endl*/;
    /*cout << "Final number of cells at G1 = " << getNumG1() << endl*/;
    /*cout << "Final number of cells at S = " << getNumS() << endl*/;
    /*cout << "Final number of cells at G2 = " << getNumG2() << endl*/;
    /*cout << "Final number of cells at M = " << getNumM() << endl*/;
    /*cout << "Final number of cells at G0 = " << getNumG0() << endl*/;
    /*cout << "Final number of living cells = "
         << getNumFib() << endl*/;
    /*cout << "Final number of tumor cells = " << getNumTum() << endl*/;
    /*cout << "Final tumor density: " << OUT_TUM_DENS << "%" << endl*/;
    /*cout << (PAR_INIT_TUM_DENS - OUT_TUM_DENS) / PAR_INIT_TUM_DENS
            * 100 << "% of initial tumor cells killed" << endl*/;
    /*cout << "Final number of vessels = " << getNumVes() << endl*/;
    /*cout << "Final vascular density: " << OUT_VES_DENS << "%" << endl*/;
    /*cout << "Final number of dead cells = "<< getNumDead() << endl*/;
    /*cout << "---------------------------------------------" << endl*/;
    return 0;
}


int Tissue::updateModel(const double currentTime,
                        const double DT){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->updateModel(currentTime, DT);
    }

    ST_PREV_TUM_DENS = ST_TUM_DENS;
    ST_TUM_DENS      = double(getNumTum()) / double(m_numComp) * 100.0;
    ST_INT_TUM_DENS += 0.5 * DT * (ST_PREV_TUM_DENS + ST_TUM_DENS);

    if(m_treatment){
        double tumSurv;
        tumSurv = ST_TUM_DENS / PAR_INIT_TUM_DENS;
        if(tumSurv < 0.5 && m_doseNeeded[0] < 0.0){
            m_doseNeeded[0] = ((Cell*)m_comp->at(0))->getAccDose();
            m_timeNeeded[0] = currentTime;
            printNeededDoseAndTime(0);
        }
        if(tumSurv < 0.2 && m_doseNeeded[1] < 0.0){
            m_doseNeeded[1] = ((Cell*)m_comp->at(0))->getAccDose();
            m_timeNeeded[1] = currentTime;
            printNeededDoseAndTime(1);
        }
        if(tumSurv < 0.1 && m_doseNeeded[2] < 0.0){
            m_doseNeeded[2] = ((Cell*)m_comp->at(0))->getAccDose();
            m_timeNeeded[2] = currentTime;
            printNeededDoseAndTime(2);
        }
        if(tumSurv < 0.05 && m_doseNeeded[3] < 0.0){
            m_doseNeeded[3] = ((Cell*)m_comp->at(0))->getAccDose();
            m_timeNeeded[3] = currentTime;
            printNeededDoseAndTime(3);
        }
        if(tumSurv < 0.01 && m_doseNeeded[4] < 0.0){
            m_doseNeeded[4] = ((Cell*)m_comp->at(0))->getAccDose();
            m_timeNeeded[4] = currentTime;
            printNeededDoseAndTime(4);
        }
        if(tumSurv < 0.001 && m_doseNeeded[5] < 0.0){
            m_doseNeeded[5] = ((Cell*)m_comp->at(0))->getAccDose();
            m_timeNeeded[5] = currentTime;
            printNeededDoseAndTime(5);
        }
    }
    return 0;
}


int Tissue::getNumFib() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getFib()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumDead() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getDead()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumG0() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getG0()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumG1() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getG1()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumG2() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getG2()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumM() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getM()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumNormVes() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getNormVes()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumS() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getS()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumTum() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getTum()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumTumVes() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getTumVes()){
            count++;
        }
    }
    return count;
}


int Tissue::getNumVes() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((Cell *)m_comp->at(k))->getVes()){
            count++;
        }
    }
    return count;
}


Treatment *Tissue::getTreatment() const{
    return m_treatment;
}


void Tissue::printNeededDoseAndTime(int sel) const{
    vector<string> perc = {"50", "80", "90", "95", "99", "99.9"};

    /*/*cout << "Total dose needed to kill " << perc[sel] <<
            "% of tumor cells = " << m_doseNeeded[sel] <<
            " Gy" << endl*/;
    /*cout << "Time needed to kill " << perc[sel] <<
            "% of tumor cells = " << m_timeNeeded[sel] <<
            " h" << endl*/;
}

