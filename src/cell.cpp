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

Cell::Cell(Model *const parent) : Model(10, 17, 2, 36, 0){
    ST_FIB     = 1.0;
    ST_TUM     = 0.0;
    ST_VES     = 0.0;
    ST_TUM_VES = 0.0;
    ST_HYP_NEC = 0.0;
    ST_MIT_CAT = 0.0;
    ST_APOP    = 0.0;
    ST_DEAD    = 0.0;

    ST_TIMER = 0.0; //h

    PAR_TUM_GROWTH = 1.0;
    PAR_DOUB_TIME  = 1008.0; //h

    PAR_LIM_G1S = 0.55 * PAR_DOUB_TIME;
    PAR_LIM_SG2 = 0.75 * PAR_DOUB_TIME;
    PAR_LIM_G2M = 0.9  * PAR_DOUB_TIME;

    ST_G1 = 0.0;
    ST_S  = 0.0;
    ST_G2 = 0.0;
    ST_M  = 0.0;
    ST_G0 = 0.0;

    PAR_RES           = 1.0;
    PAR_FIB_DOUB_TIME = 234.0;

    PAR_ANG        = 1.0;
    PAR_ANG_TIME   = 5040.0; //h
    PAR_VEGF       = 0.0;
    PAR_VEGF_THRES = 32.5;

    PAR_ALPHA_FIB      = 0.0; //Gy^-1
    PAR_ALPHA_G1       = 0.158; //Gy^-1
    PAR_ALPHA_S        = 0.113; //Gy^-1
    PAR_ALPHA_G2       = 0.169; //Gy^-1
    PAR_ALPHA_M        = 0.189; //Gy^-1
    PAR_ALPHA_G0       = 0.189; //Gy^-1
    PAR_ALPHA_NORM_VES = 0.0; //Gy^-1
    PAR_ALPHA_TUM_VES  = 0.0; //Gy^-1

    PAR_ALPHA = PAR_ALPHA_FIB;

    PAR_BETA_FIB      = 0.0; //Gy^-2
    PAR_BETA_G1       = 0.051; //Gy^-2
    PAR_BETA_S        = 0.037; //Gy^-2
    PAR_BETA_G2       = 0.055; //Gy^-2
    PAR_BETA_M        = 0.061; //Gy^-2
    PAR_ALPHA_G0      = 0.061; //Gy^-2
    PAR_BETA_NORM_VES = 0.0; //Gy^-2
    PAR_BETA_TUM_VES  = 0.0; //Gy^-2

    PAR_BETA = PAR_BETA_FIB;

    ST_DAM = 0.0;

	PAR_DOSE_THRES  = 6.0;
    PAR_ARREST_TIME = 18.0;
    ST_ARREST = 0.0;

    PAR_M   = 3.0; //adim.
    PAR_K   = 3.0; //mmHg
    PAR_PO2 = 3.0; //mmHg
    PAR_HYP_NEC_THRES = 1.0; //mmHg

    PAR_ACC_DOSE = 0.0; //Gy

    m_parent = parent;
    m_treatment = 0;
    m_edge = new vector<Cell *>((unsigned int)0, 0);
}


Cell::Cell(const int i, const int j, const int l,
           const double tumGrowth, const double doubTime,
           vector <double> cycDur, const double res,
           const double fibDoubTime, const double ang,
           const double angTime, const double vegfThres,
           vector<double> alpha, vector<double> beta,
           const double doseThres, const double arrestTime,
           const double hypNecThres, Model *const parent) :
    Model(10, 17, 2, 36, 0){
    m_i = i;
    m_j = j;
    m_l = l;

    ST_FIB      = 1.0;
    ST_TUM      = 0.0;
    ST_NORM_VES = 0.0;
    ST_TUM_VES  = 0.0;
    ST_VES      = 0.0;
    ST_HYP_NEC  = 0.0;
    ST_MIT_CAT  = 0.0;
    ST_APOP     = 0.0;
    ST_DEAD     = 0.0;

    ST_TIMER = 0.0; //h

    PAR_TUM_GROWTH = tumGrowth;
    PAR_DOUB_TIME  = doubTime; //h

    PAR_LIM_G1S = cycDur.at(0) * PAR_DOUB_TIME;
    PAR_LIM_SG2 = PAR_LIM_G1S + cycDur.at(1) * PAR_DOUB_TIME;
    PAR_LIM_G2M = PAR_LIM_SG2 + cycDur.at(2) * PAR_DOUB_TIME;

    ST_G1 = 0.0;
    ST_S  = 0.0;
    ST_G2 = 0.0;
    ST_M  = 0.0;
    ST_G0 = 0.0;

    PAR_RES            = res;
    PAR_FIB_DOUB_TIME  = fibDoubTime;

    PAR_ANG        = ang;
    PAR_ANG_TIME   = angTime; //h
    PAR_VEGF       = 0.0;
    PAR_VEGF_THRES = vegfThres;

    PAR_ALPHA_FIB      = alpha.at(0); //Gy^-1
    PAR_ALPHA_G1       = alpha.at(1); //Gy^-1
    PAR_ALPHA_S        = alpha.at(2); //Gy^-1
    PAR_ALPHA_G2       = alpha.at(3); //Gy^-1
    PAR_ALPHA_M        = alpha.at(4); //Gy^-1
    PAR_ALPHA_G0       = alpha.at(5); //Gy^-1
    PAR_ALPHA_NORM_VES = alpha.at(6); //Gy^-1
    PAR_ALPHA_TUM_VES  = alpha.at(7); //Gy^-1

    PAR_ALPHA = PAR_ALPHA_FIB;

    PAR_BETA_FIB      = beta.at(0); //Gy^-2
    PAR_BETA_G1       = beta.at(1); //Gy^-2
    PAR_BETA_S        = beta.at(2); //Gy^-2
    PAR_BETA_G2       = beta.at(3); //Gy^-2
    PAR_BETA_M        = beta.at(4); //Gy^-2
    PAR_BETA_G0       = beta.at(5); //Gy^-2
    PAR_BETA_NORM_VES = beta.at(6); //Gy^-2
    PAR_BETA_TUM_VES  = beta.at(7); //Gy^-2

    PAR_BETA = PAR_BETA_FIB;

    ST_DAM = 0.0;

	PAR_DOSE_THRES  = doseThres;
    PAR_ARREST_TIME = arrestTime;
    ST_ARREST = 0.0;

    PAR_M   = 3.0; //adim.
    PAR_K   = 3.0; //mmHg
    PAR_PO2 = 0.0; //mmHg
    PAR_HYP_NEC_THRES = hypNecThres; //mmHg

    PAR_ACC_DOSE = 0; //Gy

    m_parent = parent;
    m_treatment = ((Tissue *)m_parent)->getTreatment();
    m_edge = new vector<Cell *>((unsigned int)0, 0);
}


Cell::~Cell(){
    delete m_edge;
}


int Cell::calcModelOut(){
    OUT_STATE = ST_FIB + 2 * (ST_TUM && !ST_DAM) +
            3 * (ST_TUM && ST_DAM) + 4 * ST_NORM_VES +
            5 * ST_TUM_VES + 6 * ST_HYP_NEC +
            7 * ST_MIT_CAT + 8 * ST_APOP;
    OUT_CYCLE = ST_G1 + 2 * ST_S + 3 * ST_G2 +
            4 * ST_M + 5 * ST_G0;
    return 0;
}


int Cell::initModel(){
    ST_FIB     = !IN_TUM && !IN_NORM_VES && !IN_TUM_VES &&
            !IN_HYP_NEC && !IN_APOP && !IN_MIT_CAT;
    ST_TUM     = IN_TUM && !IN_HYP_NEC && !IN_APOP &&
            !IN_MIT_CAT;
    ST_NORM_VES = IN_NORM_VES && !IN_HYP_NEC && !IN_APOP &&
            !IN_MIT_CAT;
    ST_TUM_VES  = IN_TUM_VES && !IN_HYP_NEC && !IN_APOP &&
            !IN_MIT_CAT;
    ST_VES      = ST_NORM_VES || ST_TUM_VES;
    ST_HYP_NEC  = IN_HYP_NEC && !ST_FIB && !ST_TUM &&
            !ST_VES;
    ST_MIT_CAT  = IN_MIT_CAT && !ST_FIB && !ST_TUM &&
            !ST_VES;
    ST_APOP     = IN_APOP && !ST_FIB && !ST_TUM &&
            !ST_VES;
    ST_DEAD     = ST_HYP_NEC || ST_APOP || ST_MIT_CAT;

    setInFib(0.0);
    setInTum(0.0);
    setInNormVes(0.0);
    setInTumVes(0.0);
    setInHypNec(0.0);
    setInMitCat(0.0);
    setInApop(0.0);

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

    if(ST_FIB){
        PAR_ALPHA = PAR_ALPHA_FIB;
        PAR_BETA  = PAR_BETA_FIB;
    }
    if(ST_TUM){
        if(ST_G1){
            PAR_ALPHA = PAR_ALPHA_G1;
            PAR_BETA  = PAR_BETA_G1;
        }
        if(ST_S){
            PAR_ALPHA = PAR_ALPHA_S;
            PAR_BETA  = PAR_BETA_S;
        }
        if(ST_G2){
            PAR_ALPHA = PAR_ALPHA_G2;
            PAR_BETA  = PAR_BETA_G2;
        }
        if(ST_M){
            PAR_ALPHA = PAR_ALPHA_M;
            PAR_BETA  = PAR_BETA_M;
        }
        if(ST_G0){
            PAR_ALPHA = PAR_ALPHA_G0;
            PAR_BETA  = PAR_BETA_G0;
        }
    }
    if(ST_NORM_VES){
        PAR_ALPHA = PAR_ALPHA_NORM_VES;
        PAR_BETA  = PAR_BETA_NORM_VES;
    }
    if(ST_TUM_VES){
        PAR_ALPHA = PAR_ALPHA_TUM_VES;
        PAR_BETA  = PAR_BETA_TUM_VES;
    }
    if(ST_DEAD){
        PAR_ALPHA = 0.0;
        PAR_BETA  = 0.0;
    }

    PAR_PO2  = IN_PO2;
    PAR_VEGF = IN_VEGF;
    return 0;
}


int Cell::terminateModel(){
    return 0;
}


int Cell::updateModel(const double currentTime,
                      const double DT){
    PAR_PO2  = IN_PO2;
    PAR_VEGF = IN_VEGF;

    if(!ST_DEAD){
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
                PAR_ACC_DOSE += m_treatment->getFraction();
                if(!ST_DEAD){
                    calcRespToIrr();
                }
            }
        }
    }

    ST_FIB      = (ST_FIB || IN_FIB) && !IN_TUM && !IN_NORM_VES &&
            !IN_TUM_VES && !IN_HYP_NEC && !IN_MIT_CAT && !IN_APOP;
    ST_TUM      = (ST_TUM || IN_TUM) && !IN_NORM_VES && !IN_TUM_VES &&
            !IN_HYP_NEC && !IN_MIT_CAT && !IN_APOP;
    ST_NORM_VES = (ST_NORM_VES || IN_NORM_VES) && !IN_HYP_NEC &&
            !IN_MIT_CAT && !IN_APOP;
    ST_TUM_VES  = (ST_TUM_VES || IN_TUM_VES) && !IN_HYP_NEC && !IN_MIT_CAT &&
            !IN_APOP;
    ST_VES      = ST_NORM_VES || ST_TUM_VES;
    ST_MIT_CAT  = (ST_MIT_CAT || IN_MIT_CAT) && !IN_HYP_NEC && !ST_FIB &&
            !ST_TUM && !ST_VES;
    ST_APOP     = (ST_APOP || IN_APOP) && !IN_HYP_NEC && !ST_FIB && !ST_TUM &&
            !ST_VES;
    ST_HYP_NEC  = (ST_HYP_NEC || IN_HYP_NEC) && !ST_FIB &&
            !ST_TUM && !ST_VES;
    ST_DEAD     = ST_HYP_NEC || ST_APOP || ST_MIT_CAT;

    setInFib(0.0);
    setInTum(0.0);
    setInNormVes(0.0);
    setInTumVes(0.0);
    setInHypNec(0.0);
    setInMitCat(0.0);
    setInApop(0.0);

    if(!ST_TUM){
        ST_G1 = 0.0;
        ST_S  = 0.0;
        ST_G2 = 0.0;
        ST_M  = 0.0;
        ST_G0 = 0.0;
    }

    if(ST_FIB){
        PAR_ALPHA = PAR_ALPHA_FIB;
        PAR_BETA  = PAR_BETA_FIB;
    }
    else if(ST_TUM){
        if(ST_G1){
            PAR_ALPHA = PAR_ALPHA_G1;
            PAR_BETA  = PAR_BETA_G1;
        }
        else if(ST_S){
            PAR_ALPHA = PAR_ALPHA_S;
            PAR_BETA  = PAR_BETA_S;
        }
        else if(ST_G2){
            PAR_ALPHA = PAR_ALPHA_G2;
            PAR_BETA  = PAR_BETA_G2;
        }
        else if(ST_M){
            PAR_ALPHA = PAR_ALPHA_M;
            PAR_BETA  = PAR_BETA_M;
        }
        else if(ST_G0){
            PAR_ALPHA = PAR_ALPHA_G0;
            PAR_BETA  = PAR_BETA_G0;
        }
    }
    else if(ST_NORM_VES){
        PAR_ALPHA = PAR_ALPHA_NORM_VES;
        PAR_BETA  = PAR_BETA_NORM_VES;
    }
    else if(ST_TUM_VES){
        PAR_ALPHA = PAR_ALPHA_TUM_VES;
        PAR_BETA  = PAR_BETA_TUM_VES;
    }
    else if(ST_DEAD){
        PAR_ALPHA = 0.0;
        PAR_BETA  = 0.0;
    }
    return 0;
}


void Cell::addToEdge(Cell *const cell){
    m_edge->push_back(cell);
}


void Cell::calcFibProlif(double DT){
    if(ST_TIMER >= PAR_FIB_DOUB_TIME){
        Cell *newFib(0);
        newFib = searchSpaceForFib();
        if(newFib){
            ST_TIMER = 0.0;
            newFib->setInFib(1.0);
            newFib->ST_TIMER = 0.0;
        }
    }

    else{
        ST_TIMER += DT;
    }
}


double Cell::calcG1SF() const{
    //double PAR_C(2.0), PAR_PO2_INFL(2.0);
    //return exp(-exp(-PAR_C * (PAR_PO2 - PAR_PO2_INFL)));
    return 1.0;
}


double Cell::calcG2MF() const{
    //double PAR_C(2.0), PAR_PO2_INFL(2.0);
    //return exp(-exp(-PAR_C * (PAR_PO2 - PAR_PO2_INFL)));
    return 1.0;
}


void Cell::calcHypNec(){
    if(PAR_PO2 < PAR_HYP_NEC_THRES){
        setInHypNec(1.0);
    }
}


void Cell::calcNormVesProlif(double DT){
    if(ST_TIMER >= PAR_ANG_TIME){
        if(PAR_VEGF >= PAR_VEGF_THRES){
            Cell *newVes(0);
            newVes = searchSpaceForTumVes();
            if(newVes){
                ST_TIMER = 0.0;
                newVes->setInTumVes(1.0);
                newVes->ST_TIMER = 0.0;
            }
        }
    }

    else{
        ST_TIMER += DT;
    }
}


double Cell::calcOER() const{
    return (PAR_M * PAR_PO2 + PAR_K) / (PAR_PO2 + PAR_K); //mmHg
}


void Cell::calcRespToIrr(){
    ST_ARREST = PAR_ARREST_TIME;

    double p(double(rand()) / double(RAND_MAX));
    if(calcSF() < p){
        if(m_treatment->getFraction() < PAR_DOSE_THRES){
            if(PAR_TUM_GROWTH){
                ST_DAM = 1.0;
            }
            else{
                setInApop(1.0);
            }
        }
        else{
            setInApop(1.0);
        }
    }
}


double Cell::calcSF() const{
    double fraction, OER, SF;

    fraction = m_treatment->getFraction();
    OER = calcOER();
    SF = exp(-PAR_ALPHA / PAR_M * fraction * OER -
             PAR_BETA / (PAR_M * PAR_M) * fraction * fraction *
             OER * OER);
    return SF;
}


void Cell::calcTumGrowth(double DT){
    if(ST_G1 && ST_TIMER >= PAR_LIM_G1S){
        if(ST_ARREST > 0.0){
            ST_ARREST -= DT;
        }
        else{
            ST_TIMER += DT;
            ST_G1 = 0.0;
            ST_S  = 1.0;
        }
    }

    else if(ST_S && ST_TIMER >= PAR_LIM_SG2){
        ST_TIMER += DT;
        ST_S  = 0.0;
        ST_G2 = 1.0;
    }

    else if(ST_G2 && ST_TIMER >= PAR_LIM_G2M){
        if(ST_ARREST > 0.0){
            ST_ARREST -= DT;
        }
        else{
            ST_TIMER += DT;
            ST_G2 = 0.0;
            ST_M  = 1.0;
        }
    }

    else if(ST_TIMER >= PAR_DOUB_TIME){
        Cell *newTumCell(0);
        newTumCell = searchSpaceForTum();
        if(newTumCell){
            ST_TIMER = 0.0;
            if(ST_DAM){
                ST_DAM = 0.0;
                setInMitCat(1.0);
            }
            else{
                ST_M  = 0.0;
                ST_G0 = 0.0;
                ST_G1 = 1.0;
                newTumCell->setInTum(1.0);
                newTumCell->ST_TIMER = 0.0;
                newTumCell->ST_G1 = 1.0;
                newTumCell->ST_ARREST = ST_ARREST;
            }
        }
        else{
            ST_M  = 0.0;
            ST_G0 = 1.0;
        }
    }

    else{
        ST_TIMER += DT;
    }
}


void Cell::calcTumVesProlif(double DT){
    if(ST_TIMER >= PAR_ANG_TIME){
        if(PAR_VEGF >= PAR_VEGF_THRES){
            Cell *newVes(0);
            newVes = searchSpaceForTumVes();
            if(newVes){
                ST_TIMER = 0.0;
                newVes->setInTumVes(1.0);
                newVes->ST_TIMER = 0.0;
            }
        }
    }
    else{
        ST_TIMER += DT;
    }
}


bool Cell::getFib() const{
    return ST_FIB;
}


double Cell::getAccDose() const{
    return PAR_ACC_DOSE;
}


bool Cell::getDead() const{
    return ST_DEAD;
}


double Cell::getDoubTime() const{
    return PAR_DOUB_TIME;
}


vector<Cell *> *Cell::getEdge() const{
    return m_edge;
}


bool Cell::getG0() const{
    return ST_G0;
}


bool Cell::getG1() const{
    return ST_G1;
}


bool Cell::getG2() const{
    return ST_G2;
}


bool Cell::getHypNec() const{
    return ST_HYP_NEC;
}


bool Cell::getM() const{
    return ST_M;
}


bool Cell::getMitCat() const{
    return ST_MIT_CAT;
}


bool Cell::getNormVes() const{
    return ST_NORM_VES;
}


int Cell::getOutState() const{
    return int(OUT_STATE);
}


double Cell::getR() const{
    return m_r;
}


bool Cell::getS() const{
    return ST_S;
}


bool Cell::getTum() const{
    return ST_TUM;
}


bool Cell::getTumVes() const{
    return ST_TUM_VES;
}


bool Cell::getTumDam() const{
    return ST_TUM && ST_DAM;
}


bool Cell::getTumNotDam() const{
    return ST_TUM && !ST_DAM;
}


bool Cell::getVes() const{
    return ST_VES;
}


Cell *Cell::searchInitSpaceForTum() const{
    int edgeSize, m;

    edgeSize = m_edge->size();
    m = rand() % edgeSize;
    for(int n(0); n < edgeSize; n++){
        if(!m_edge->at(m)->IN_TUM &&
                !m_edge->at(m)->IN_NORM_VES &&
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


Cell *Cell::searchSpaceForFib() const{
    int edgeSize, m;

    edgeSize = m_edge->size();
    m = rand() % edgeSize;
    for(int n(0); n < edgeSize; n++){
        if((m_edge->at(m)->getOutState() == 6 ||
            m_edge->at(m)->getOutState() == 7 ||
            m_edge->at(m)->getOutState() == 8) &&
                !m_edge->at(m)->IN_FIB &&
                !m_edge->at(m)->IN_TUM &&
                !m_edge->at(m)->IN_NORM_VES &&
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


Cell *Cell::searchSpaceForTum() const{
    int edgeSize, m;

    edgeSize = m_edge->size();
    m = rand() % edgeSize;
    for(int n(0); n < edgeSize; n++){
        if((m_edge->at(m)->getOutState() == 1 ||
            m_edge->at(m)->getOutState() == 7 ||
            m_edge->at(m)->getOutState() == 8) &&
                !m_edge->at(m)->IN_FIB &&
                !m_edge->at(m)->IN_TUM &&
                !m_edge->at(m)->IN_NORM_VES &&
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


Cell *Cell::searchSpaceForTumVes() const{
    int edgeSize;
    double m(0.0);
    Cell *newTumVes(0), *tempCell(0);

    edgeSize = m_edge->size();

    for(int n(0); n < edgeSize; n++){
        tempCell = m_edge->at(n);
        if(tempCell->getOutState() != 4 &&
                tempCell->getOutState() != 5 &&
                tempCell->getOutState() != 6 &&
                !tempCell->IN_FIB &&
                !tempCell->IN_TUM &&
                !tempCell->IN_NORM_VES &&
                !tempCell->IN_TUM_VES){
            if(tempCell->PAR_VEGF > m){
                m = tempCell->PAR_VEGF;
                newTumVes = tempCell;
            }
        }
    }
    return newTumVes;
}


Cell *Cell::searchSpaceForVes() const{
    return 0;
}


void Cell::setInFib(const double input){
    IN_FIB = input;
}


void Cell::setInApop(const double input){
    IN_APOP = input;
}


void Cell::setInHypNec(const double input){
    IN_HYP_NEC = input;
}


void Cell::setInMitCat(const double input){
    IN_MIT_CAT = input;
}


void Cell::setInNormVes(const double input){
    IN_NORM_VES = input;
}


void Cell::setInPO2(const double input){
    IN_PO2 = input;
}


void Cell::setInTimer(const double input){
    IN_TIMER = input;
}


void Cell::setInTum(const double input){
    IN_TUM = input;
}


void Cell::setInTumVes(const double input){
    IN_TUM_VES = input;
}


void Cell::setInVegf(const double input){
    IN_VEGF = input;
}


void Cell::setR(const Cell *origCell){
    m_r = sqrt((m_i - origCell->m_i) * (m_i - origCell->m_i) +
               (m_j - origCell->m_j) * (m_j - origCell->m_j) +
               (m_l - origCell->m_l) * (m_l - origCell->m_l));
}


bool compRedge(Cell *a, Cell *b){
    return (a->getR() < b->getR());
}


