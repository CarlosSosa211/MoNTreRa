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
                     const double Dvegf, const double DO2, const double Vmax,
                     const double Km, const double pO2NormVes,
                     const double pO2TumVes, const double hypThres,
                     const double VmaxVegf, const double KmVegf,
                     const double hypVegf) :
    Model(0, 0, 5, 2, nrow * ncol * nlayer){
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

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new OxyCell(Vmax, Km, pO2NormVes, pO2TumVes,
                                    hypThres, VmaxVegf, KmVegf, hypVegf,
                                    this);
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
                     const double Dvegf, const double DO2, const double Vmax,
                     const double Km, const double pO2NormVes,
                     const double pO2TumVes, const double hypThres,
                     const double  VmaxVegf, const double KmVegf,
                     const double hypVegf) :
    Model(0, 0, 5, 2, nrow * ncol * nlayer){
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

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new OxyCell(Vmax, Km, pO2NormVes, pO2TumVes,
                                    hypThres, VmaxVegf, KmVegf, hypVegf,
                                    this);
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
    for (int i(0); i < m_numComp; i++){
        ((OxyCell *)(m_comp->at(i)))->initModel();
    }
    return 1;
}


int OxyTissue::calcModelOut(){
    OUT_HYP_DENS = double(getNumHyp()) / double(m_numComp) * 100.0;

    vector<double> pO2, vegf;
    for(int k(0); k < m_numComp; k++){
        (m_comp->at(k))->calcModelOut();
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

    return 0;
}


int OxyTissue::updateModel(double currentTime, const double DT){
    if(m_nlayer == 1){
        //First row, first column
        m_map[0][0][0]->
                setInDiffO2(PAR_DO2 * (m_map[0][0][1]->getOutPO2() +
                m_map[0][1][0]->getOutPO2() -
                2.0 * m_map[0][0][0]->getOutPO2()));

        m_map[0][0][0]->
                setInDiffVEGF(PAR_DVEGF * (m_map[0][0][1]->getOutVEGF() +
                m_map[0][1][0]->getOutVEGF() -
                2.0 * m_map[0][0][0]->getOutVEGF()));

        //First row, last column
        m_map[0][0][m_ncol - 1]->
                setInDiffO2(PAR_DO2 * (m_map[0][0][m_ncol - 2]->getOutPO2() +
                m_map[0][1][m_ncol - 1]->getOutPO2() -
                2.0 * m_map[0][0][m_ncol - 1]->getOutPO2()));

        m_map[0][0][m_ncol - 1]->
                setInDiffVEGF(PAR_DVEGF * (m_map[0][0][m_ncol - 2]->getOutVEGF() +
                m_map[0][1][m_ncol - 1]->getOutVEGF() -
                2.0 * m_map[0][0][m_ncol - 1]->getOutVEGF()));

        //Last row, first column
        m_map[0][m_nrow - 1][0]->
                setInDiffO2(PAR_DO2 * (m_map[0][m_nrow - 2][0]->getOutPO2() +
                m_map[0][m_nrow - 1][1]->getOutPO2() -
                2.0 * m_map[0][m_nrow - 1][0]->getOutPO2()));

        m_map[0][m_nrow - 1][0]->
                setInDiffVEGF(PAR_DVEGF * (m_map[0][m_nrow - 2][0]->getOutVEGF() +
                m_map[0][m_nrow - 1][1]->getOutVEGF() -
                2.0 * m_map[0][m_nrow - 1][0]->getOutVEGF()));

        //Last row, last column
        m_map[0][m_nrow - 1][m_ncol - 1]->
                setInDiffO2(PAR_DO2 * (m_map[0][m_nrow - 2][m_ncol - 1]->getOutPO2() +
                m_map[0][m_nrow - 1][m_ncol - 2]->getOutPO2() -
                2.0 * m_map[0][m_nrow - 1][m_ncol - 1]->getOutPO2()));

        m_map[0][m_nrow - 1][m_ncol - 1]->
                setInDiffVEGF(PAR_DVEGF * (m_map[0][m_nrow - 2][m_ncol - 1]->getOutVEGF() +
                m_map[0][m_nrow - 1][m_ncol - 2]->getOutVEGF() -
                2.0 * m_map[0][m_nrow - 1][m_ncol - 1]->getOutVEGF()));

        //First row, middle columns
        for(int j(1); j < m_ncol - 1; j++){
            m_map[0][0][j]->
                    setInDiffO2(PAR_DO2 * (m_map[0][0][j - 1]->getOutPO2() +
                    m_map[0][0][j + 1]->getOutPO2() +
                    m_map[0][1][j]->getOutPO2() -
                    3.0 * m_map[0][0][j]->getOutPO2()));

            m_map[0][0][j]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[0][0][j - 1]->getOutVEGF() +
                    m_map[0][0][j + 1]->getOutVEGF() +
                    m_map[0][1][j]->getOutVEGF() -
                    3.0 * m_map[0][0][j]->getOutVEGF()));
        }

        //Last row, middle columns
        for(int j(1); j < m_ncol - 1; j++){
            m_map[0][m_nrow - 1][j]->
                    setInDiffO2(PAR_DO2 * (m_map[0][m_nrow - 1][j - 1]->getOutPO2() +
                    m_map[0][m_nrow - 1][j + 1]->getOutPO2() +
                    m_map[0][m_nrow - 2][j]->getOutPO2() -
                    3.0 * m_map[0][m_nrow - 1][j]->getOutPO2()));

            m_map[0][m_nrow - 1][j]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[0][m_nrow - 1][j - 1]->getOutVEGF() +
                    m_map[0][m_nrow - 1][j + 1]->getOutVEGF() +
                    m_map[0][m_nrow - 2][j]->getOutVEGF() -
                    3.0 * m_map[0][m_nrow - 1][j]->getOutVEGF()));
        }

        //First column, middle rows
        for(int i(1); i < m_nrow - 1; i++){
            m_map[0][i][0]->
                    setInDiffO2(PAR_DO2 * (m_map[0][i - 1][0]->getOutPO2() +
                    m_map[0][i][1]->getOutPO2() +
                    m_map[0][i + 1][0]->getOutPO2() -
                    3.0 * m_map[0][i][0]->getOutPO2()));

            m_map[0][i][0]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[0][i - 1][0]->getOutVEGF() +
                    m_map[0][i][1]->getOutVEGF() +
                    m_map[0][i + 1][0]->getOutVEGF() -
                    3.0 * m_map[0][i][0]->getOutVEGF()));
        }

        //Last column, middle rows
        for(int i(1); i < m_nrow - 1; i++){
            m_map[0][i][m_ncol - 1]->
                    setInDiffO2(PAR_DO2 * (m_map[0][i - 1][m_ncol - 1]->getOutPO2() +
                    m_map[0][i][m_ncol - 2]->getOutPO2() +
                    m_map[0][i + 1][m_ncol - 1]->getOutPO2() -
                    3.0 * m_map[0][i][m_ncol - 1]->getOutPO2()));

            m_map[0][i][m_ncol - 1]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[0][i - 1][m_ncol - 1]->getOutVEGF() +
                    m_map[0][i][m_ncol - 2]->getOutVEGF() +
                    m_map[0][i + 1][m_ncol - 1]->getOutVEGF() -
                    3.0 * m_map[0][i][m_ncol - 1]->getOutVEGF()));
        }

        //Middle rows, middle columns
        for(int i(1); i < m_nrow - 1; i++){
            for(int j(1); j < m_ncol - 1; j++){
                m_map[0][i][j]->
                        setInDiffO2(PAR_DO2 * (m_map[0][i - 1][j]->getOutPO2() +
                        m_map[0][i][j - 1]->getOutPO2() +
                        m_map[0][i][j + 1]->getOutPO2() +
                        m_map[0][i + 1][j]->getOutPO2() -
                        4.0 * m_map[0][i][j]->getOutPO2()));

                m_map[0][i][j]->
                        setInDiffVEGF(PAR_DVEGF * (m_map[0][i - 1][j]->getOutVEGF() +
                        m_map[0][i][j - 1]->getOutVEGF() +
                        m_map[0][i][j + 1]->getOutVEGF() +
                        m_map[0][i + 1][j]->getOutVEGF() -
                        4.0 * m_map[0][i][j]->getOutVEGF()));
            }
        }
    }

    else{
        //First layer, first row, first column
        m_map[0][0][0]->
                setInDiffO2(PAR_DO2 * (m_map[0][0][1]->getOutPO2() +
                m_map[0][1][0]->getOutPO2() +
                m_map[1][0][0]->getOutPO2() -
                3.0 * m_map[0][0][0]->getOutPO2()));

        m_map[0][0][0]->
                setInDiffVEGF(PAR_DVEGF * (m_map[0][0][1]->getOutVEGF() +
                m_map[0][1][0]->getOutVEGF() +
                m_map[1][0][0]->getOutVEGF() -
                3.0 * m_map[0][0][0]->getOutVEGF()));

        //First layer, first row, last column
        m_map[0][0][m_ncol - 1]->
                setInDiffO2(PAR_DO2 * (m_map[0][0][m_ncol - 2]->getOutPO2() +
                m_map[0][1][m_ncol - 1]->getOutPO2() +
                m_map[1][0][m_ncol - 1]->getOutPO2() -
                3.0 * m_map[0][0][m_ncol - 1]->getOutPO2()));

        m_map[0][0][m_ncol - 1]->
                setInDiffVEGF(PAR_DVEGF * (m_map[0][0][m_ncol - 2]->getOutVEGF() +
                m_map[0][1][m_ncol - 1]->getOutVEGF() +
                m_map[1][0][m_ncol - 1]->getOutVEGF() -
                3.0 * m_map[0][0][m_ncol - 1]->getOutVEGF()));

        //First layer, Last row, first column
        m_map[0][m_nrow - 1][0]->
                setInDiffO2(PAR_DO2 * (m_map[0][m_nrow - 2][0]->getOutPO2() +
                m_map[0][m_nrow - 1][1]->getOutPO2() +
                m_map[1][m_nrow - 1][0]->getOutPO2() -
                3.0 * m_map[0][m_nrow - 1][0]->getOutPO2()));

        m_map[0][m_nrow - 1][0]->
                setInDiffVEGF(PAR_DVEGF * (m_map[0][m_nrow - 2][0]->getOutVEGF() +
                m_map[0][m_nrow - 1][1]->getOutVEGF() +
                m_map[1][m_nrow - 1][0]->getOutVEGF() -
                3.0 * m_map[0][m_nrow - 1][0]->getOutVEGF()));

        //First layer, last row, last column
        m_map[0][m_nrow - 1][m_ncol - 1]->
                setInDiffO2(PAR_DO2 * (m_map[0][m_nrow - 2][m_ncol - 1]->getOutPO2() +
                m_map[0][m_nrow - 1][m_ncol - 2]->getOutPO2() +
                m_map[1][m_nrow - 1][m_ncol - 1]->getOutPO2() -
                3.0 * m_map[0][m_nrow - 1][m_ncol - 1]->getOutPO2()));

        m_map[0][m_nrow - 1][m_ncol - 1]->
                setInDiffVEGF(PAR_DVEGF * (m_map[0][m_nrow - 2][m_ncol - 1]->getOutVEGF() +
                m_map[0][m_nrow - 1][m_ncol - 2]->getOutVEGF() +
                m_map[1][m_nrow - 1][m_ncol - 1]->getOutVEGF() -
                3.0 * m_map[0][m_nrow - 1][m_ncol - 1]->getOutVEGF()));

        //Last layer, first row, first column
        m_map[m_nlayer - 1][0][0]->
                setInDiffO2(PAR_DO2 * (m_map[m_nlayer - 1][0][1]->getOutPO2() +
                m_map[m_nlayer - 1][1][0]->getOutPO2() +
                m_map[m_nlayer - 2][0][0]->getOutPO2() -
                3.0 * m_map[m_nlayer - 1][0][0]->getOutPO2()));

        m_map[m_nlayer - 1][0][0]->
                setInDiffVEGF(PAR_DVEGF * (m_map[m_nlayer - 1][0][1]->getOutVEGF() +
                m_map[m_nlayer - 1][1][0]->getOutVEGF() +
                m_map[m_nlayer - 2][0][0]->getOutVEGF() -
                3.0 * m_map[m_nlayer - 1][0][0]->getOutVEGF()));

        //Last layer, first row, last column
        m_map[m_nlayer - 1][0][m_ncol - 1]->
                setInDiffO2(PAR_DO2 * (m_map[m_nlayer - 1][0][m_ncol - 2]->getOutPO2() +
                m_map[m_nlayer - 1][1][m_ncol - 1]->getOutPO2() +
                m_map[m_nlayer - 2][0][m_ncol - 1]->getOutPO2() -
                3.0 * m_map[m_nlayer - 1][0][m_ncol - 1]->getOutPO2()));

        m_map[m_nlayer - 1][0][m_ncol - 1]->
                setInDiffVEGF(PAR_DVEGF * (m_map[m_nlayer - 1][0][m_ncol - 2]->getOutVEGF() +
                m_map[m_nlayer - 1][1][m_ncol - 1]->getOutVEGF() +
                m_map[m_nlayer - 2][0][m_ncol - 1]->getOutVEGF() -
                3.0 * m_map[m_nlayer - 1][0][m_ncol - 1]->getOutVEGF()));

        //Last layer, Last row, first column
        m_map[m_nlayer - 1][m_nrow - 1][0]->
                setInDiffO2(PAR_DO2 * (m_map[m_nlayer - 1][m_nrow - 2][0]->getOutPO2() +
                m_map[m_nlayer - 1][m_nrow - 1][1]->getOutPO2() +
                m_map[m_nlayer - 2][m_nrow - 1][0]->getOutPO2() -
                3.0 * m_map[m_nlayer - 1][m_nrow - 1][0]->getOutPO2()));

        m_map[m_nlayer - 1][m_nrow - 1][0]->
                setInDiffVEGF(PAR_DVEGF * (m_map[m_nlayer - 1][m_nrow - 2][0]->getOutVEGF() +
                m_map[m_nlayer - 1][m_nrow - 1][1]->getOutVEGF() +
                m_map[m_nlayer - 2][m_nrow - 1][0]->getOutVEGF() -
                3.0 * m_map[m_nlayer - 1][m_nrow - 1][0]->getOutVEGF()));

        //Last layer, last row, last column
        m_map[m_nlayer - 1][m_nrow - 1][m_ncol - 1]->
                setInDiffO2(PAR_DO2 * (m_map[m_nlayer - 1][m_nrow - 2][m_ncol - 1]->getOutPO2() +
                m_map[m_nlayer - 1][m_nrow - 1][m_ncol - 2]->getOutPO2() +
                m_map[m_nrow - 2][m_nrow - 1][m_ncol - 1]->getOutPO2() -
                3.0 * m_map[m_nrow - 1][m_nrow - 1][m_ncol - 1]->getOutPO2()));

        m_map[m_nlayer - 1][m_nrow - 1][m_ncol - 1]->
                setInDiffVEGF(PAR_DVEGF * (m_map[m_nlayer - 1][m_nrow - 2][m_ncol - 1]->getOutVEGF() +
                m_map[m_nlayer - 1][m_nrow - 1][m_ncol - 2]->getOutVEGF() +
                m_map[m_nlayer - 2][m_nrow - 1][m_ncol - 1]->getOutVEGF() -
                3.0 * m_map[m_nlayer - 1][m_nrow - 1][m_ncol - 1]->getOutVEGF()));

        //First layer, first row, middle columns
        for(int j(1); j < m_ncol - 1; j++){
            m_map[0][0][j]->
                    setInDiffO2(PAR_DO2 * (m_map[0][0][j - 1]->getOutPO2() +
                    m_map[0][0][j + 1]->getOutPO2() +
                    m_map[0][1][j]->getOutPO2() +
                    m_map[1][0][j]->getOutPO2() -
                    4.0 * m_map[0][0][j]->getOutPO2()));

            m_map[0][0][j]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[0][0][j - 1]->getOutVEGF() +
                    m_map[0][0][j + 1]->getOutVEGF() +
                    m_map[0][1][j]->getOutVEGF() +
                    m_map[1][0][j]->getOutVEGF() -
                    4.0 * m_map[0][0][j]->getOutVEGF()));
        }

        //First layer, last row, middle columns
        for(int j(1); j < m_ncol - 1; j++){
            m_map[0][m_nrow - 1][j]->
                    setInDiffO2(PAR_DO2 * (m_map[0][m_nrow - 1][j - 1]->getOutPO2() +
                    m_map[0][m_nrow - 1][j + 1]->getOutPO2() +
                    m_map[0][m_nrow - 2][j]->getOutPO2() +
                    m_map[1][m_nrow - 1][j]->getOutPO2() -
                    4.0 * m_map[0][m_nrow - 1][j]->getOutPO2()));

            m_map[0][m_nrow - 1][j]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[0][m_nrow - 1][j - 1]->getOutVEGF() +
                    m_map[0][m_nrow - 1][j + 1]->getOutVEGF() +
                    m_map[0][m_nrow - 2][j]->getOutVEGF() +
                    m_map[1][m_nrow - 1][j]->getOutVEGF() -
                    4.0 * m_map[0][m_nrow - 1][j]->getOutVEGF()));
        }

        //First layer, first column, middle rows
        for(int i(1); i < m_nrow - 1; i++){
            m_map[0][i][0]->
                    setInDiffO2(PAR_DO2 * (m_map[0][i - 1][0]->getOutPO2() +
                    m_map[0][i][1]->getOutPO2() +
                    m_map[0][i + 1][0]->getOutPO2() +
                    m_map[1][i][0]->getOutPO2() -
                    4.0 * m_map[0][i][0]->getOutPO2()));

            m_map[0][i][0]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[0][i - 1][0]->getOutVEGF() +
                    m_map[0][i][1]->getOutVEGF() +
                    m_map[0][i + 1][0]->getOutVEGF() +
                    m_map[1][i][0]->getOutVEGF() -
                    4.0 * m_map[0][i][0]->getOutVEGF()));
        }

        //First layer, last column, middle rows
        for(int i(1); i < m_nrow - 1; i++){
            m_map[0][i][m_ncol - 1]->
                    setInDiffO2(PAR_DO2 * (m_map[0][i - 1][m_ncol - 1]->getOutPO2() +
                    m_map[0][i][m_ncol - 2]->getOutPO2() +
                    m_map[0][i + 1][m_ncol - 1]->getOutPO2() +
                    m_map[1][i][m_ncol - 1]->getOutPO2() -
                    4.0 * m_map[0][i][m_ncol - 1]->getOutPO2()));

            m_map[0][i][m_ncol - 1]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[0][i - 1][m_ncol - 1]->getOutVEGF() +
                    m_map[0][i][m_ncol - 2]->getOutVEGF() +
                    m_map[0][i + 1][m_ncol - 1]->getOutVEGF() +
                    m_map[1][i][m_ncol - 1]->getOutVEGF() -
                    4.0 * m_map[0][i][m_ncol - 1]->getOutVEGF()));
        }

        //Last layer, first row, middle columns
        for(int j(1); j < m_ncol - 1; j++){
            m_map[m_nlayer - 1][0][j]->
                    setInDiffO2(PAR_DO2 * (m_map[m_nlayer - 1][0][j - 1]->getOutPO2() +
                    m_map[m_nlayer - 1][0][j + 1]->getOutPO2() +
                    m_map[m_nlayer - 1][1][j]->getOutPO2() +
                    m_map[m_nlayer - 2][0][j]->getOutPO2() -
                    4.0 * m_map[m_nlayer - 1][0][j]->getOutPO2()));

            m_map[m_nlayer - 1][0][j]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[m_nlayer - 1][0][j - 1]->getOutVEGF() +
                    m_map[m_nlayer - 1][0][j + 1]->getOutVEGF() +
                    m_map[m_nlayer - 1][1][j]->getOutVEGF() +
                    m_map[m_nlayer - 2][0][j]->getOutVEGF() -
                    4.0 * m_map[m_nlayer - 1][0][j]->getOutVEGF()));
        }

        //Last layer, last row, middle columns
        for(int j(1); j < m_ncol - 1; j++){
            m_map[m_nlayer - 1][m_nrow - 1][j]->
                    setInDiffO2(PAR_DO2 *  (m_map[m_nlayer - 1][m_nrow - 1][j - 1]->getOutPO2() +
                    m_map[m_nlayer - 1][m_nrow - 1][j + 1]->getOutPO2() +
                    m_map[m_nlayer - 1][m_nrow - 2][j]->getOutPO2() +
                    m_map[m_nlayer - 2][m_nrow - 1][j]->getOutPO2() -
                    4.0 * m_map[m_nlayer - 1][m_nrow - 1][j]->getOutPO2()));

            m_map[m_nlayer - 1][m_nrow - 1][j]->
                    setInDiffVEGF(PAR_DVEGF *  (m_map[m_nlayer - 1][m_nrow - 1][j - 1]->getOutVEGF() +
                    m_map[m_nlayer - 1][m_nrow - 1][j + 1]->getOutVEGF() +
                    m_map[m_nlayer - 1][m_nrow - 2][j]->getOutVEGF() +
                    m_map[m_nlayer - 2][m_nrow - 1][j]->getOutVEGF() -
                    4.0 * m_map[m_nlayer - 1][m_nrow - 1][j]->getOutVEGF()));
        }

        //Last layer, first column, middle rows
        for(int i(1); i < m_nrow - 1; i++){
            m_map[m_nlayer - 1][i][0]->
                    setInDiffO2(PAR_DO2 * (m_map[m_nlayer - 1][i - 1][0]->getOutPO2() +
                    m_map[m_nlayer - 1][i][1]->getOutPO2() +
                    m_map[m_nlayer - 1][i + 1][0]->getOutPO2() +
                    m_map[m_nlayer - 2][i][0]->getOutPO2() -
                    4.0 * m_map[m_nlayer - 1][i][0]->getOutPO2()));

            m_map[m_nlayer - 1][i][0]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[m_nlayer - 1][i - 1][0]->getOutVEGF() +
                    m_map[m_nlayer - 1][i][1]->getOutVEGF() +
                    m_map[m_nlayer - 1][i + 1][0]->getOutVEGF() +
                    m_map[m_nlayer - 2][i][0]->getOutVEGF() -
                    4.0 * m_map[m_nlayer - 1][i][0]->getOutVEGF()));
        }

        //Last layer, last column, middle rows
        for(int i(1); i < m_nrow - 1; i++){
            m_map[m_nlayer - 1][i][m_ncol - 1]->
                    setInDiffO2(PAR_DO2 * (m_map[m_nlayer - 1][i - 1][m_ncol - 1]->getOutPO2() +
                    m_map[m_nlayer - 1][i][m_ncol - 2]->getOutPO2() +
                    m_map[m_nlayer - 1][i + 1][m_ncol - 1]->getOutPO2() +
                    m_map[m_nlayer - 2][i][m_ncol - 1]->getOutPO2() -
                    4.0 * m_map[m_nlayer - 1][i][m_ncol - 1]->getOutPO2()));

            m_map[m_nlayer - 1][i][m_ncol - 1]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[m_nlayer - 1][i - 1][m_ncol - 1]->getOutVEGF() +
                    m_map[m_nlayer - 1][i][m_ncol - 2]->getOutVEGF() +
                    m_map[m_nlayer - 1][i + 1][m_ncol - 1]->getOutVEGF() +
                    m_map[m_nlayer - 2][i][m_ncol - 1]->getOutVEGF() -
                    4.0 * m_map[m_nlayer - 1][i][m_ncol - 1]->getOutVEGF()));
        }

        //Middle layers, first row, first column
        for(int l(1); l < m_nlayer- 1; l++){
            m_map[l][0][0]->
                    setInDiffO2(PAR_DO2 * (m_map[l][0][1]->getOutPO2() +
                    m_map[l][1][0]->getOutPO2() +
                    m_map[l + 1][0][0]->getOutPO2() +
                    m_map[l - 1][0][0]->getOutPO2() -
                    4.0 * m_map[l][0][0]->getOutPO2()));

            m_map[l][0][0]->
                    setInDiffVEGF(PAR_DVEGF *  (m_map[l][0][1]->getOutVEGF() +
                    m_map[l][1][0]->getOutVEGF() +
                    m_map[l + 1][0][0]->getOutVEGF() +
                    m_map[l - 1][0][0]->getOutVEGF() -
                    4.0 * m_map[l][0][0]->getOutVEGF()));
        }

        //Middle layers, first row, last column
        for(int l(1); l < m_nlayer - 1; l++){
            m_map[l][0][m_ncol - 1]->
                    setInDiffO2(PAR_DO2 * (m_map[l][0][m_ncol - 2]->getOutPO2() +
                    m_map[l][1][m_ncol - 1]->getOutPO2() +
                    m_map[l + 1][0][m_ncol - 1]->getOutPO2() +
                    m_map[l - 1][0][m_ncol - 1]->getOutPO2() -
                    4.0 * m_map[l][0][m_ncol - 1]->getOutPO2()));

            m_map[l][0][m_ncol - 1]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[l][0][m_ncol - 2]->getOutVEGF() +
                    m_map[l][1][m_ncol - 1]->getOutVEGF() +
                    m_map[l + 1][0][m_ncol - 1]->getOutVEGF() +
                    m_map[l - 1][0][m_ncol - 1]->getOutVEGF() -
                    4.0 * m_map[l][0][m_ncol - 1]->getOutVEGF()));
        }

        //Middle layers, last row, first column
        for(int l(1); l < m_nlayer - 1; l++){
            m_map[l][m_nrow - 1][0]->
                    setInDiffO2(PAR_DO2 * (m_map[l][m_nrow - 1][1]->getOutPO2() +
                    m_map[l][m_nrow - 2][0]->getOutPO2() +
                    m_map[l + 1][m_nrow - 1][0]->getOutPO2() +
                    m_map[l - 1][m_nrow - 1][0]->getOutPO2() -
                    4.0 * m_map[l][m_nrow - 1][0]->getOutPO2()));

            m_map[l][m_nrow - 1][0]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[l][m_nrow - 1][1]->getOutVEGF() +
                    m_map[l][m_nrow - 2][0]->getOutVEGF() +
                    m_map[l + 1][m_nrow - 1][0]->getOutVEGF() +
                    m_map[l - 1][m_nrow - 1][0]->getOutVEGF() -
                    4.0 * m_map[l][m_nrow - 1][0]->getOutVEGF()));
        }

        //Middle layers, last row, last column
        for(int l(1); l < m_nlayer - 1; l++){
            m_map[l][m_nrow - 1][m_ncol - 1]->
                    setInDiffO2(PAR_DO2 * (m_map[l][m_nrow - 1][m_ncol - 2]->getOutPO2() +
                    m_map[l][m_nrow - 2][m_ncol - 1]->getOutPO2() +
                    m_map[l + 1][m_nrow - 1][m_ncol - 1]->getOutPO2() +
                    m_map[l - 1][m_nrow - 1][m_ncol - 1]->getOutPO2() -
                    4.0 * m_map[l][m_nrow - 1][m_ncol - 1]->getOutPO2()));

            m_map[l][m_nrow - 1][m_ncol - 1]->
                    setInDiffVEGF(PAR_DVEGF * (m_map[l][m_nrow - 1][m_ncol - 2]->getOutVEGF() +
                    m_map[l][m_nrow - 2][m_ncol - 1]->getOutVEGF() +
                    m_map[l + 1][m_nrow - 1][m_ncol - 1]->getOutVEGF() +
                    m_map[l - 1][m_nrow - 1][m_ncol - 1]->getOutVEGF() -
                    4.0 * m_map[l][m_nrow - 1][m_ncol - 1]->getOutVEGF()));
        }

        //First layer, middle rows, middle columns
        for(int i(1); i < m_nrow - 1; i++){
            for(int j(1); j < m_ncol - 1; j++){
                m_map[0][i][j]->
                        setInDiffO2(PAR_DO2 * (m_map[0][i - 1][j]->getOutPO2() +
                        m_map[0][i][j - 1]->getOutPO2() +
                        m_map[0][i][j + 1]->getOutPO2() +
                        m_map[0][i + 1][j]->getOutPO2() +
                        m_map[1][i][j]->getOutPO2() -
                        5.0 * m_map[0][i][j]->getOutPO2()));

                m_map[0][i][j]->
                        setInDiffVEGF(PAR_DVEGF * (m_map[0][i - 1][j]->getOutVEGF() +
                        m_map[0][i][j - 1]->getOutVEGF() +
                        m_map[0][i][j + 1]->getOutVEGF() +
                        m_map[0][i + 1][j]->getOutVEGF() +
                        m_map[1][i][j]->getOutVEGF() -
                        5.0 * m_map[0][i][j]->getOutVEGF()));
            }
        }

        //Last layer, middle rows, middle columns
        for(int i(1); i < m_nrow - 1; i++){
            for(int j(1); j < m_ncol - 1; j++){
                m_map[m_nlayer - 1][i][j]->
                        setInDiffO2(PAR_DO2 * (m_map[m_nlayer - 1][i - 1][j]->getOutPO2() +
                        m_map[m_nlayer - 1][i][j - 1]->getOutPO2() +
                        m_map[m_nlayer - 1][i][j + 1]->getOutPO2() +
                        m_map[m_nlayer - 1][i + 1][j]->getOutPO2() +
                        m_map[m_nlayer - 2][i][j]->getOutPO2() -
                        5.0 * m_map[m_nlayer - 1][i][j]->getOutPO2()));

                m_map[m_nlayer - 1][i][j]->
                        setInDiffVEGF(PAR_DVEGF * (m_map[m_nlayer - 1][i - 1][j]->getOutVEGF() +
                        m_map[m_nlayer - 1][i][j - 1]->getOutVEGF() +
                        m_map[m_nlayer - 1][i][j + 1]->getOutVEGF() +
                        m_map[m_nlayer - 1][i + 1][j]->getOutVEGF() +
                        m_map[m_nlayer - 2][i][j]->getOutVEGF() -
                        5.0 * m_map[m_nlayer - 1][i][j]->getOutVEGF()));
            }
        }

        //Middle layers, middle rows, first column
        for(int l(1); l < m_nlayer - 1; l++){
            for(int i(1); i < m_nrow - 1; i++){
                m_map[l][i][0]->
                        setInDiffO2(PAR_DO2 * (m_map[l + 1][i][0]->getOutPO2() +
                        m_map[l - 1][i][0]->getOutPO2() +
                        m_map[l][i + 1][0]->getOutPO2() +
                        m_map[l][i - 1][0]->getOutPO2() +
                        m_map[l][i][1]->getOutPO2() -
                        5.0 * m_map[l][i][0]->getOutPO2()));

                m_map[l][i][0]->
                        setInDiffVEGF(PAR_DVEGF * (m_map[l + 1][i][0]->getOutVEGF() +
                        m_map[l - 1][i][0]->getOutVEGF() +
                        m_map[l][i + 1][0]->getOutVEGF() +
                        m_map[l][i - 1][0]->getOutVEGF() +
                        m_map[l][i][1]->getOutVEGF() -
                        5.0 * m_map[l][i][0]->getOutVEGF()));
            }
        }

        //Middle layers, middle rows, last column
        for(int l(1); l < m_nlayer - 1; l++){
            for(int i(1); i < m_nrow - 1; i++){
                m_map[l][i][m_ncol - 1]->
                        setInDiffO2(PAR_DO2 * (m_map[l + 1][i][m_ncol - 1]->getOutPO2() +
                        m_map[l - 1][i][m_ncol - 1]->getOutPO2() +
                        m_map[l][i + 1][m_ncol - 1]->getOutPO2() +
                        m_map[l][i - 1][m_ncol - 1]->getOutPO2() +
                        m_map[l][i][m_ncol - 2]->getOutPO2() -
                        5.0 * m_map[l][i][m_ncol - 1]->getOutPO2()));

                m_map[l][i][m_ncol - 1]->
                        setInDiffVEGF(PAR_DVEGF * (m_map[l + 1][i][m_ncol - 1]->getOutVEGF() +
                        m_map[l - 1][i][m_ncol - 1]->getOutVEGF() +
                        m_map[l][i + 1][m_ncol - 1]->getOutVEGF() +
                        m_map[l][i - 1][m_ncol - 1]->getOutVEGF() +
                        m_map[l][i][m_ncol - 2]->getOutVEGF() -
                        5.0 * m_map[l][i][m_ncol - 1]->getOutVEGF()));
            }
        }

        //Middle layers, first rows, middle columns
        for(int l(1); l < m_nlayer - 1; l++){
            for(int j(1); j < m_ncol - 1; j++){
                m_map[l][0][j]->
                        setInDiffO2(PAR_DO2 * (m_map[l + 1][0][j]->getOutPO2() +
                        m_map[l - 1][0][j]->getOutPO2() +
                        m_map[l][0][j + 1]->getOutPO2() +
                        m_map[l][0][j - 1]->getOutPO2() +
                        m_map[l][1][j]->getOutPO2() -
                        5.0 * m_map[l][0][j]->getOutPO2()));

                m_map[l][0][j]->
                        setInDiffVEGF(PAR_DVEGF * (m_map[l + 1][0][j]->getOutVEGF() +
                        m_map[l - 1][0][j]->getOutVEGF() +
                        m_map[l][0][j + 1]->getOutVEGF() +
                        m_map[l][0][j - 1]->getOutVEGF() +
                        m_map[l][1][j]->getOutVEGF() -
                        5.0 * m_map[l][0][j]->getOutVEGF()));
            }
        }

        //Middle layers, last rows, middle columns
        for(int l(1); l < m_nlayer - 1; l++){
            for(int j(1); j < m_ncol - 1; j++){
                m_map[l][m_nrow - 1][j]->
                        setInDiffO2(PAR_DO2 * (m_map[l + 1][m_nrow - 1][j]->getOutPO2() +
                        m_map[l - 1][m_nrow - 1][j]->getOutPO2() +
                        m_map[l][m_nrow - 1][j + 1]->getOutPO2() +
                        m_map[l][m_nrow - 1][j - 1]->getOutPO2() +
                        m_map[l][m_nrow - 2][j]->getOutPO2() -
                        5.0 * m_map[l][m_nrow - 1][j]->getOutPO2()));

                m_map[l][m_nrow - 1][j]->
                        setInDiffVEGF(PAR_DVEGF * (m_map[l + 1][m_nrow - 1][j]->getOutVEGF() +
                        m_map[l - 1][m_nrow - 1][j]->getOutVEGF() +
                        m_map[l][m_nrow - 1][j + 1]->getOutVEGF() +
                        m_map[l][m_nrow - 1][j - 1]->getOutVEGF() +
                        m_map[l][m_nrow - 2][j]->getOutVEGF() -
                        5.0 * m_map[l][m_nrow - 1][j]->getOutVEGF()));
            }
        }

        //Middle layers, middle rows, middle columns
        for(int l(1); l < m_nlayer - 1; l++){
            for(int i(1); i < m_nrow - 1; i++){
                for(int j(1); j < m_ncol - 1; j++){
                    m_map[l][i][j]->
                            setInDiffO2(PAR_DO2 * (m_map[l + 1][i][j]->getOutPO2() +
                                        m_map[l - 1][i][j]->getOutPO2() +
                            m_map[l][i + 1][j]->getOutPO2() +
                            m_map[l][i - 1][j]->getOutPO2() +
                            m_map[l][i][j + 1]->getOutPO2() +
                            m_map[l][i][j - 1]->getOutPO2() -
                            6.0 * m_map[l][i][j]->getOutPO2()));

                    m_map[l][i][j]->
                            setInDiffVEGF(PAR_DVEGF * (m_map[l + 1][i][j]->getOutVEGF() +
                                          m_map[l - 1][i][j]->getOutVEGF() +
                            m_map[l][i + 1][j]->getOutVEGF() +
                            m_map[l][i - 1][j]->getOutVEGF() +
                            m_map[l][i][j + 1]->getOutVEGF() +
                            m_map[l][i][j - 1]->getOutVEGF() -
                            6.0 * m_map[l][i][j]->getOutVEGF()));
                }
            }
        }
    }

    for(int i(0); i < m_numComp; i++){
        (m_comp->at(i))->updateModel(currentTime, DT);
    }


    /*Not sure that this is working. Members in if may be equal.*/
    bool stable(true);
    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j++){
                if(fabs(m_map[l][i][j]->getOutPO2() -
                        m_map[l][i][j]->getPO2()) >
                        1e-2 * m_map[l][i][j]->getOutPO2() ||
                        fabs(m_map[l][i][j]->getOutVEGF() -
                             m_map[l][i][j]->getVEGF()) >
                        1e-2 * m_map[l][i][j]->getOutVEGF()){
                    stable = false;
                }
            }
        }
    }
    if(stable){
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
