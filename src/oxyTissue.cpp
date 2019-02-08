/**
 * @file oxyTissue.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17.0
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

#include "oxyTissue.hpp"

using namespace std;

OxyTissue::OxyTissue(const int nrow, const int ncol, const int nlayer,
                     const double cellSize, const string nFInVes,
                     const double ang, const double Dvegf, const double VmaxVegf,
                     const double KmVegf, const double hypVegf,
                     const double DO2, const double VmaxO2, const double KmO2,
                     const double pO2NormVes, const double pO2TumVes,
                     const double hypThres) :
    Model(0, 4, 9, 3, nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    if(m_nlayer == 1){
        PAR_DO2   = DO2 / (cellSize * cellSize);
        PAR_DVEGF = Dvegf / (cellSize * cellSize);
    }

    else{
        PAR_DO2   = DO2 / (cellSize * cellSize * cellSize);
        PAR_DVEGF = Dvegf / (cellSize * cellSize * cellSize);
    }

    PAR_OXY_ANG = ang;

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new OxyCell(ang, VmaxVegf, KmVegf,
                                    hypVegf, VmaxO2, KmO2,
                                    pO2NormVes, pO2TumVes,
                                    hypThres, this);
        m_numOut += (m_comp->at(k))->getNumOut();
    }


    m_map = new OxyCell ***[m_nlayer];
    for(int l(0); l < m_nlayer; l++){
        m_map[l] = new OxyCell **[m_nrow];
        for(int i(0); i < m_nrow; i++){
            m_map[l][i] = new OxyCell *[m_ncol];
        }
    }

    int k(0);
    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j ++){
                m_map[l][i][j] = ((OxyCell *)m_comp->at(k));
                k++;
            }
        }
    }

    int ii, jj, ll;
    const int tii[6] = {0, 0, 0, -1, 1, 0};
    const int tjj[6] = {0, -1, 1, 0, 0, 0};
    const int tll[6] = {-1, 0, 0, 0, 0, 1};

    int incol, lnrowNcol, m, mm;

    for(int l(0); l < m_nlayer; l++){
        lnrowNcol = l * m_nrow * m_ncol;
        for(int i(0); i < m_nrow; i++){
            incol = i * m_ncol;
            for(int j(0); j < m_ncol; j++){
                m = lnrowNcol + incol + j;
                for(int n(0); n < 6; n++){
                    ii = tii[n];
                    jj = tjj[n];
                    ll = tll[n];
                    mm = lnrowNcol + ll * m_nrow * m_ncol + incol + ii * m_ncol + j + jj;
                    if(l + ll >= 0 && l + ll < m_nlayer && i + ii >= 0 && i + ii < m_nrow &&
                            j + jj >= 0 && j + jj < m_ncol){
                        ((OxyCell *)m_comp->at(m))->addToEdge((OxyCell *)m_comp->at(mm));
                    }
                }
            }
        }
    }

    ifstream fInVes(nFInVes.c_str());
    double inputVes;

    for(int k(0); k < m_numComp; k++){
        if(fInVes >> inputVes){
            ((OxyCell *)m_comp->at(k))->setInNormVes(inputVes);
        }
        else{
            cout << "Insufficient data in vessel file" << endl;
            break;
        }
    }
    fInVes.close();
}


OxyTissue::OxyTissue(const int nrow, const int ncol, const int nlayer,
                     const double cellSize, const vector<bool> &inVes,
                     const double ang, const double Dvegf, const double VmaxVegf,
                     const double KmVegf, const double hypVegf,
                     const double DO2, const double VmaxO2, const double KmO2,
                     const double pO2NormVes, const double pO2TumVes,
                     const double hypThres) :
    Model(0, 4, 9, 3, nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    if(m_nlayer == 1){
        PAR_DO2   = DO2 / (cellSize * cellSize);
        PAR_DVEGF = Dvegf / (cellSize * cellSize);
    }

    else{
        PAR_DO2   = DO2 / (cellSize * cellSize * cellSize);
        PAR_DVEGF = Dvegf / (cellSize * cellSize * cellSize);
    }

    PAR_OXY_ANG = ang;

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new OxyCell(ang, VmaxVegf, KmVegf,
                                    hypVegf, VmaxO2, KmO2,
                                    pO2NormVes, pO2TumVes,
                                    hypThres, this);
        m_numOut += (m_comp->at(k))->getNumOut();
        ((OxyCell *)m_comp->at(k))->setInNormVes(inVes.at(k));
    }

    m_map = new OxyCell ***[m_nlayer];
    for(int l(0); l < m_nlayer; l++){
        m_map[l] = new OxyCell **[m_nrow];
        for(int i(0); i < m_nrow; i++){
            m_map[l][i] = new OxyCell *[m_ncol];
        }
    }

    int k(0);
    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j ++){
                m_map[l][i][j] = ((OxyCell *)m_comp->at(k));
                k++;
            }
        }
    }

    int ii, jj, ll;
    int incol, lnrowNcol, m, mm;
    const int tii[6] = {0, 0, 0, -1, 1, 0};
    const int tjj[6] = {0, -1, 1, 0, 0, 0};
    const int tll[6] = {-1, 0, 0, 0, 0, 1};

    for(int l(0); l < m_nlayer; l++){
        lnrowNcol = l * m_nrow * m_ncol;
        for(int i(0); i < m_nrow; i++){
            incol = i * m_ncol;
            for(int j(0); j < m_ncol; j++){
                m = lnrowNcol + incol + j;
                for(int n(0); n < 6; n++){
                    ii = tii[n];
                    jj = tjj[n];
                    ll = tll[n];
                    mm = lnrowNcol + ll * m_nrow * m_ncol + incol + ii * m_ncol + j + jj;
                    if(l + ll >= 0 && l + ll < m_nlayer && i + ii >= 0 && i + ii < m_nrow &&
                            j + jj >= 0 && j + jj < m_ncol){
                        ((OxyCell *)m_comp->at(m))->addToEdge(((OxyCell *)m_comp->at(mm)));
                    }
                }
            }
        }
    }
}


OxyTissue::~OxyTissue(){
    for(int k(0); k < m_numComp; k++){
        delete m_comp->at(k);
    }

    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            delete [] m_map[l][i];
        }
        delete [] m_map[l];
    }
    delete [] m_map;
}


int OxyTissue::initModel(){
    ST_OXYSTABLE  = 0.0;
    ST_VEGFSTABLE = 0.0;
    ST_TIME_TO_OXYSTABLE  = 0.0;
    ST_TIME_TO_VEGFSTABLE = 0.0;

    for (int k(0); k < m_numComp; k++){
        ((OxyCell *)(m_comp->at(k)))->initModel();
    }
    return 0;
}


int OxyTissue::calcModelOut(){
    OUT_HYP_DENS = double(getNumHyp()) / double(m_numComp) * 100.0;

    vector<double> pO2, vegf;
    for(int k(0); k < m_numComp; k++){
        m_comp->at(k)->calcModelOut();
        pO2.push_back(m_comp->at(k)->getOut()->at(0));
        vegf.push_back(m_comp->at(k)->getOut()->at(1));
    }
    int npO2 (pO2.size() / 2), nvegf (vegf.size() / 2);
    nth_element(pO2.begin(), pO2.begin() + npO2, pO2.end());
    nth_element(vegf.begin(), vegf.begin() + nvegf, vegf.end());
    OUT_PO2_MED  = pO2.at(npO2);
    OUT_VEGF_MED = vegf.at(nvegf);

    OUT_PO2_MEAN  = accumulate(pO2.begin(), pO2.end(), 0.0) / pO2.size();
    OUT_VEGF_MEAN = accumulate(vegf.begin(), vegf.end(), 0.0) / vegf.size();

    OUT_OXYSTABLE  = ST_OXYSTABLE;
    OUT_VEGFSTABLE = ST_VEGFSTABLE;
    OUT_TIME_TO_OXYSTABLE  = ST_TIME_TO_OXYSTABLE;
    OUT_TIME_TO_VEGFSTABLE = ST_TIME_TO_VEGFSTABLE;

    return 0;
}


int OxyTissue::updateModel(double currentTime, const double DT){
    int edgeSize;
    double diffO2, diffVegf;
    std::vector<OxyCell *> *edge;

    if(!ST_OXYSTABLE){
        for(int l(0); l < m_nlayer; l++){
            for(int i(0); i < m_nrow; i++){
                for(int j(0); j < m_ncol; j++){
                    edge = m_map[l][i][j]->getEdge();
                    edgeSize = edge->size();
                    diffO2   = 0.0;
                    for(int n(0); n < edgeSize; n++){
                        diffO2   += edge->at(n)->getOutPO2();
                    }
                    diffO2   -= edgeSize * m_map[l][i][j]->getOutPO2();
                    m_map[l][i][j]->setInDiffO2(PAR_DO2 * diffO2);
                }
            }
        }
    }

    if(!ST_VEGFSTABLE){
        for(int l(0); l < m_nlayer; l++){
            for(int i(0); i < m_nrow; i++){
                for(int j(0); j < m_ncol; j++){
                    edge = m_map[l][i][j]->getEdge();
                    edgeSize = edge->size();
                    diffVegf = 0.0;
                    for(int n(0); n < edgeSize; n++){
                        diffVegf += edge->at(n)->getOutVEGF();
                    }
                    diffVegf -= edgeSize * m_map[l][i][j]->getOutVEGF();
                    m_map[l][i][j]->setInDiffVEGF(PAR_DVEGF * diffVegf);
                }
            }
        }
    }

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k)->updateModel(currentTime, DT);
    }

    ST_OXYSTABLE  = isOxyStable();
    ST_VEGFSTABLE = isVegfStable();

    if(ST_OXYSTABLE && !ST_TIME_TO_OXYSTABLE){
        ST_TIME_TO_OXYSTABLE = currentTime;
    }

    if(ST_VEGFSTABLE && !ST_TIME_TO_VEGFSTABLE){
        ST_TIME_TO_VEGFSTABLE = currentTime;
    }

    if(ST_OXYSTABLE && ST_OXYVEGF){
        ST_OXYSTABLE  = 0.0;
        ST_VEGFSTABLE = 0.0;
        return 1;
    }
    else{
        return 0;
    }
}


int OxyTissue::terminateModel(){
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->terminateModel();
    }
    return 0;
}


int OxyTissue::getNumHyp() const{
    int count(0);
    for(int k(0); k < m_numComp; k++){
        if(((OxyCell *)m_comp->at(k))->getHyp()){
            count++;
        }
    }
    return count;
}

bool OxyTissue::isOxyStable() const{
    int k(0);
    bool stable(true);

    while(stable && k < m_numComp){
        stable = ((OxyCell *)m_comp->at(k))->getOxyStable();
        k++;
    }
    return stable;
}


bool OxyTissue::isVegfStable() const{
    int k(0);
    bool stable(true);

    while(stable && k < m_numComp){
        stable = ((OxyCell *)m_comp->at(k))->getVegfStable();
        k++;
    }
    return stable;
}
