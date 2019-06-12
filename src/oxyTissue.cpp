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

/*------------------------------------------------------------------------------
 * Constructor of the class OxyTissue.
 *
 * Inputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue (um),
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - ang: angiogenesis,
 *  - DVegf: VEGF diffusion coefficient (um^2/ms),
 *  - VmaxVegf: maximum VEGF consumption ratio (mol/um^3ms),
 *  - KmVegf: VEGF Michaelis constant (mol/um^3),
 *  - hypVegf: VEGF concentration fixed value for hypoxic cells (mol/um^3),
 *  - oxy: integer indicating the oxygenation scenario (1, space-and-time
 *  dependent),
 *  - DO2: pO2 diffusion coefficient (um^2/ms),
 *  - VmaxO2: maximum pO2 consumption ratio (mmHg/ms),
 *  - KmO2: pO2 Michaelis constant (mmHg),
 *  - pO2NormVes: fixed pO2 value for pre-existing endothelial cells (mmHg),
 *  - pO2TumVes: fixed pO2 value for neo-created endothelial cells (mmHg),
 *  - hypThres: pO2 hypoxia threshold (mmHg).
------------------------------------------------------------------------------*/

OxyTissue::OxyTissue(const int nrow, const int ncol, const int nlayer,
                     const double cellSize, const vector<bool> &inVes,
                     const bool ang, const double DVegf, const double VmaxVegf,
                     const double KmVegf, const double hypVegf, const int oxy,
                     const double DO2, const double VmaxO2, const double KmO2,
                     const double pO2NormVes, const double pO2TumVes,
                     const double hypThres) :
    AbsOxyTissue(OXYTISSUE_NUM_IN_B, OXYTISSUE_NUM_IN_I, OXYTISSUE_NUM_IN_D,
                 OXYTISSUE_NUM_ST_B, OXYTISSUE_NUM_ST_I, OXYTISSUE_NUM_ST_D,
                 OXYTISSUE_NUM_OUT_B, OXYTISSUE_NUM_OUT_I, OXYTISSUE_NUM_OUT_D,
                 OXYTISSUE_NUM_PAR_B, OXYTISSUE_NUM_PAR_I, OXYTISSUE_NUM_PAR_D,
                 nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    if(m_nlayer == 1){
        PAR_DO2   = DO2 / (cellSize * cellSize);
        PAR_DVEGF = DVegf / (cellSize * cellSize);
    }

    else{
        PAR_DO2   = DO2 / (cellSize * cellSize * cellSize);
        PAR_DVEGF = DVegf / (cellSize * cellSize * cellSize);
    }

    PAR_OXYTISSUE_ANG = ang;
    PAR_OXYTISSUE_OXY = oxy;

    int count(0);
    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new OxyCell(ang, VmaxVegf, KmVegf, hypVegf, VmaxO2,
                                    KmO2, pO2NormVes, pO2TumVes, hypThres,
                                    this);
        ((OxyCell *)m_comp->at(k))->setInNormVes(inVes.at(k));
        if(inVes.at(k)){
            count++;
        }
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
                    mm = lnrowNcol + ll * m_nrow * m_ncol + incol + ii *
                            m_ncol + j + jj;
                    if(l + ll >= 0 && l + ll < m_nlayer && i + ii >= 0 && i +
                            ii < m_nrow && j + jj >= 0 && j + jj < m_ncol){
                        ((OxyCell *)m_comp->at(m))->
                                addToEdge(((OxyCell *)m_comp->at(mm)));
                    }
                }
            }
        }
    }
}


/*------------------------------------------------------------------------------
 * Destructor of the class OxyTissue.
------------------------------------------------------------------------------*/

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


/*------------------------------------------------------------------------------
 * Redefinition of the Model initModel method.
------------------------------------------------------------------------------*/

int OxyTissue::initModel(){
    ST_OXY_STABLE  = false;
    ST_VEGF_STABLE = false;
    ST_TIME_TO_OXY_STABLE  = 0.0;
    ST_TIME_TO_VEGF_STABLE = 0.0;

    for (int k(0); k < m_numComp; k++){
        ((OxyCell *)(m_comp->at(k)))->initModel();
    }

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model calcModelOut method.
------------------------------------------------------------------------------*/

int OxyTissue::calcModelOut(){
    OUT_HYP_DENS = double(getNumHyp()) / double(m_numComp) * 100.0;

    vector<double> pO2, vegf;
    for(int k(0); k < m_numComp; k++){
        m_comp->at(k)->calcModelOut();
        pO2.push_back(m_comp->at(k)->getOutD()[0]);
        vegf.push_back(m_comp->at(k)->getOutD()[1]);
    }
    int npO2 (pO2.size() / 2), nvegf (vegf.size() / 2);
    nth_element(pO2.begin(), pO2.begin() + npO2, pO2.end());
    nth_element(vegf.begin(), vegf.begin() + nvegf, vegf.end());
    OUT_PO2_MED  = pO2.at(npO2);
    OUT_VEGF_MED = vegf.at(nvegf);

    OUT_PO2_MEAN  = accumulate(pO2.begin(), pO2.end(), 0.0) / pO2.size();
    OUT_VEGF_MEAN = accumulate(vegf.begin(), vegf.end(), 0.0) / vegf.size();

    OUT_OXYSTABLE  = ST_OXY_STABLE;
    OUT_VEGFSTABLE = ST_VEGF_STABLE;
    OUT_TIME_TO_OXYSTABLE  = ST_TIME_TO_OXY_STABLE;
    OUT_TIME_TO_VEGFSTABLE = ST_TIME_TO_VEGF_STABLE;

    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model updateModel method.
 *
 * Inputs:
 *  - currentTime: simulation current time (ms),
 *  - DT: simulation timestep (ms).
------------------------------------------------------------------------------*/

int OxyTissue::updateModel(double currentTime, const double DT){
    int edgeSize;
    double diffO2, diffVegf;
    std::vector<OxyCell *> *edge;

    //  cout << getNumVes() << endl;

    if(!ST_OXY_STABLE){
        for(int l(0); l < m_nlayer; l++){
            for(int i(0); i < m_nrow; i++){
                for(int j(0); j < m_ncol; j++){
                    edge = m_map[l][i][j]->getEdge();
                    edgeSize = edge->size();
                    diffO2   = 0.0;
                    for(int n(0); n < edgeSize; n++){
                        diffO2 += edge->at(n)->getOutPO2();
                    }
                    diffO2 -= edgeSize * m_map[l][i][j]->getOutPO2();
                    m_map[l][i][j]->setInDiffO2(PAR_DO2 * diffO2);
                }
            }
        }
    }

    if(PAR_OXYTISSUE_ANG && !ST_VEGF_STABLE){
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

    if(!ST_OXY_STABLE || !ST_VEGF_STABLE){
        for(int k(0); k < m_numComp; k++){
            m_comp->at(k)->updateModel(currentTime, DT);
        }
    }

    ST_OXY_STABLE  = isOxyStable();
    ST_VEGF_STABLE = isVegfStable();

    if(ST_OXY_STABLE && !ST_TIME_TO_OXY_STABLE){
        ST_TIME_TO_OXY_STABLE = currentTime;
    }

    if(ST_OXY_STABLE && ST_VEGF_STABLE && !ST_TIME_TO_VEGF_STABLE){
        ST_TIME_TO_VEGF_STABLE = currentTime;
    }

    if(ST_OXY_STABLE && ST_VEGF_STABLE){
        if(PAR_OXYTISSUE_OXY == 1){
            ST_OXY_STABLE  = false;
            ST_VEGF_STABLE = false;
        }
        return 1;
    }
    else{
        return 0;
    }
}


/*------------------------------------------------------------------------------
 * This function checks if the pO2 values of every cell of the tissue are
 * stable.
 *
 * Outputs:
 *  - stable: indicator of the stability of the pO2 values of the tissue
------------------------------------------------------------------------------*/

bool OxyTissue::isOxyStable() const{
    int k(0);
    bool stable(true);

    while(stable && k < m_numComp){
        stable = ((OxyCell *)m_comp->at(k))->getOxyStable();
        k++;
    }
    return stable;
}


/*------------------------------------------------------------------------------
 * This function checks if the VEGF concentration values of every cell of the
 * tissue are stable.
 *
 * Outputs:
 *  - stable: indicator of the stability of the VEGF concentration values of the
 *  tissue
------------------------------------------------------------------------------*/

bool OxyTissue::isVegfStable() const{
    int k(0);
    bool stable(true);

    while(stable && k < m_numComp){
        stable = ((OxyCell *)m_comp->at(k))->getVegfStable();
        k++;
    }
    return stable;
}
