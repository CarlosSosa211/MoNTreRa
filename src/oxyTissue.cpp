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
                     const string nFInVes, const double Dvegf,
                     const double DO2, const double Vmax, const double Km,
                     const double pO2NormVes, const double pO2TumVes,
                     const double hypThres, const double VmaxVegf,
                     const double KmVegf, const double hypVegf) :
    Model(0, 0, 5, 2, nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    PAR_DO2   = DO2;
    PAR_DVEGF = Dvegf;

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new OxyCell(Vmax, Km, pO2NormVes, pO2TumVes,
                                    hypThres, VmaxVegf, KmVegf, hypVegf,
                                    this);
        m_numOut += (m_comp->at(k))->getNumOut();
    }

    for(int l(0); l < m_nlayer; l++){
        m_map.push_back(vector<vector<Model *> >());
        for(int i(0); i < m_nrow; i++){
            m_map[l].push_back(vector<Model *>());
            for(int j(0); j < m_ncol; j ++){
                m_map[l][i].push_back(0);
            }
        }
    }

    int k(0);
    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j ++){
                m_map[l][i][j] = m_comp->at(k);
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
                     const vector<bool> &inVes, const double Dvegf,
                     const double DO2, const double Vmax, const double Km,
                     const double pO2NormVes, const double pO2TumVes,
                     const double hypThres, const double  VmaxVegf,
                     const double KmVegf, const double hypVegf) :
    Model(0, 0, 5, 2, nrow * ncol * nlayer){
    m_nrow   = nrow;
    m_ncol   = ncol;
    m_nlayer = nlayer;

    PAR_DO2   = DO2;
    PAR_DVEGF = Dvegf;

    for(int k(0); k < m_numComp; k++){
        m_comp->at(k) = new OxyCell(Vmax, Km, pO2NormVes, pO2TumVes,
                                    hypThres, VmaxVegf, KmVegf, hypVegf,
                                    this);
        m_numOut += (m_comp->at(k))->getNumOut();
        ((OxyCell *)m_comp->at(k))->setInNormVes(inVes.at(k));
    }

    for(int l(0); l < m_nlayer; l++){
        m_map.push_back(vector<vector<Model *> >());
        for(int i(0); i < m_nrow; i++){
            m_map[l].push_back(vector<Model *>());
            for(int j(0); j < m_ncol; j ++){
                m_map[l][i].push_back(0);
            }
        }
    }

    int k(0);
    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j ++){
                m_map[l][i][j] = m_comp->at(k);
                k++;
            }
        }
    }
}


OxyTissue::~OxyTissue(){
    for(int k(0); k < m_numComp; k++){
        delete m_comp->at(k);
    }
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
    double diffO2[m_nlayer][m_nrow][m_ncol];
    double diffVegf[m_nlayer][m_nrow][m_ncol];

    //First row, first column
    diffO2[0][0][0] = PAR_DO2 *
            (((OxyCell *)m_map[0][0][1])->getOutPO2() +
            ((OxyCell *)m_map[0][1][0])->getOutPO2() -
            2.0 * ((OxyCell *)m_map[0][0][0])->getOutPO2());

    diffVegf[0][0][0] = PAR_DVEGF *
            (((OxyCell *)m_map[0][0][1])->getOutVEGF() +
            ((OxyCell *)m_map[0][1][0])->getOutVEGF() -
            2.0 * ((OxyCell *)m_map[0][0][0])->getOutVEGF());

    //First row, last column
    diffO2[0][0][m_ncol - 1] = PAR_DO2 *
            (((OxyCell *)m_map[0][0][m_ncol - 2])->getOutPO2() +
            ((OxyCell *)m_map[0][1][m_ncol - 1])->getOutPO2() -
            2.0 * ((OxyCell *)m_map[0][0][m_ncol - 1])->getOutPO2());

    diffVegf[0][0][m_ncol - 1] = PAR_DVEGF *
            (((OxyCell *)m_map[0][0][m_ncol - 2])->getOutVEGF() +
            ((OxyCell *)m_map[0][1][m_ncol - 1])->getOutVEGF() -
            2.0 * ((OxyCell *)m_map[0][0][m_ncol - 1])->getOutVEGF());

    //Last row, first column
    diffO2[0][m_nrow - 1][0] = PAR_DO2 *
            (((OxyCell *)m_map[0][m_nrow - 2][0])->getOutPO2() +
            ((OxyCell *)m_map[0][m_nrow - 1][1])->getOutPO2() -
            2.0 * ((OxyCell *)m_map[0][m_nrow - 1][0])->getOutPO2());

    diffVegf[0][m_nrow - 1][0] = PAR_DVEGF *
            (((OxyCell *)m_map[0][m_nrow - 2][0])->getOutVEGF() +
            ((OxyCell *)m_map[0][m_nrow - 1][1])->getOutVEGF() -
            2.0 * ((OxyCell *)m_map[0][m_nrow - 1][0])->getOutVEGF());

    //Last row, last column
    diffO2[0][m_nrow - 1][m_ncol - 1] = PAR_DO2 *
            (((OxyCell *)m_map[0][m_nrow - 2][m_ncol - 1])->getOutPO2() +
            ((OxyCell *)m_map[0][m_nrow - 1][m_ncol - 2])->getOutPO2() -
            2.0 * ((OxyCell *)m_map[0][m_nrow - 1][m_ncol - 1])->getOutPO2());

    diffVegf[0][m_nrow - 1][m_ncol - 1] = PAR_DVEGF *
            (((OxyCell *)m_map[0][m_nrow - 2][m_ncol - 1])->getOutVEGF() +
            ((OxyCell *)m_map[0][m_nrow - 1][m_ncol - 2])->getOutVEGF() -
            2.0 * ((OxyCell *)m_map[0][m_nrow - 1][m_ncol - 1])->getOutVEGF());

    //First row, middle columns
    for(int j(1); j < m_ncol - 1; j++){
        diffO2[0][0][j] = PAR_DO2 *
                (((OxyCell *)m_map[0][0][j - 1])->getOutPO2() +
                ((OxyCell *)m_map[0][0][j + 1])->getOutPO2() +
                ((OxyCell *)m_map[0][1][j])->getOutPO2() -
                3.0 * ((OxyCell *)m_map[0][0][j])->getOutPO2());

        diffVegf[0][0][j] = PAR_DVEGF *
                (((OxyCell *)m_map[0][0][j - 1])->getOutVEGF() +
                ((OxyCell *)m_map[0][0][j + 1])->getOutVEGF() +
                ((OxyCell *)m_map[0][1][j])->getOutVEGF() -
                3.0 * ((OxyCell *)m_map[0][0][j])->getOutVEGF());
    }

    //Last row, middle columns
    for(int j(1); j < m_ncol - 1; j++){
        diffO2[0][m_nrow - 1][j] = PAR_DO2 *
                (((OxyCell *)m_map[0][m_nrow - 1][j - 1])->getOutPO2() +
                ((OxyCell *)m_map[0][m_nrow - 1][j + 1])->getOutPO2() +
                ((OxyCell *)m_map[0][m_nrow - 2][j])->getOutPO2() -
                3.0 * ((OxyCell *)m_map[0][m_nrow - 1][j])->getOutPO2());

        diffVegf[0][m_nrow - 1][j] = PAR_DVEGF *
                (((OxyCell *)m_map[0][m_nrow - 1][j - 1])->getOutVEGF() +
                ((OxyCell *)m_map[0][m_nrow - 1][j + 1])->getOutVEGF() +
                ((OxyCell *)m_map[0][m_nrow - 2][j])->getOutVEGF() -
                3.0 * ((OxyCell *)m_map[0][m_nrow - 1][j])->getOutVEGF());
    }

    //First column, middle rows
    for(int i(1); i < m_nrow - 1; i++){
        diffO2[0][i][0] = PAR_DO2 *
                (((OxyCell *)m_map[0][i - 1][0])->getOutPO2() +
                ((OxyCell *)m_map[0][i][1])->getOutPO2() +
                ((OxyCell *)m_map[0][i + 1][0])->getOutPO2() -
                3.0 * ((OxyCell *)m_map[0][i][0])->getOutPO2());

        diffVegf[0][i][0] = PAR_DVEGF *
                (((OxyCell *)m_map[0][i - 1][0])->getOutVEGF() +
                ((OxyCell *)m_map[0][i][1])->getOutVEGF() +
                ((OxyCell *)m_map[0][i + 1][0])->getOutVEGF() -
                3.0 * ((OxyCell *)m_map[0][i][0])->getOutVEGF());
    }

    //Last column, middle rows
    for(int i(1); i < m_nrow - 1; i++){
        diffO2[0][i][m_ncol - 1] = PAR_DO2 *
                (((OxyCell *)m_map[0][i - 1][m_ncol - 1])->getOutPO2() +
                ((OxyCell *)m_map[0][i][m_ncol - 2])->getOutPO2() +
                ((OxyCell *)m_map[0][i + 1][m_ncol - 1])->getOutPO2() -
                3.0 * ((OxyCell *)m_map[0][i][m_ncol - 1])->getOutPO2());

        diffVegf[0][i][m_ncol - 1] = PAR_DVEGF *
                (((OxyCell *)m_map[0][i - 1][m_ncol - 1])->getOutVEGF() +
                ((OxyCell *)m_map[0][i][m_ncol - 2])->getOutVEGF() +
                ((OxyCell *)m_map[0][i + 1][m_ncol - 1])->getOutVEGF() -
                3.0 * ((OxyCell *)m_map[0][i][m_ncol - 1])->getOutVEGF());
    }

    //Middle rows, middle columns
    for(int i(1); i < m_nrow - 1; i++){
        for(int j(1); j < m_ncol - 1; j++){
            diffO2[0][i][j] = PAR_DO2 *
                    (((OxyCell *)m_map[0][i - 1][j])->getOutPO2() +
                    ((OxyCell *)m_map[0][i][j - 1])->getOutPO2() +
                    ((OxyCell *)m_map[0][i][j + 1])->getOutPO2() +
                    ((OxyCell *)m_map[0][i + 1][j])->getOutPO2() -
                    4.0 * ((OxyCell *)m_map[0][i][j])->getOutPO2());

            diffVegf[0][i][j] = PAR_DVEGF *
                    (((OxyCell *)m_map[0][i - 1][j])->getOutVEGF() +
                    ((OxyCell *)m_map[0][i][j - 1])->getOutVEGF() +
                    ((OxyCell *)m_map[0][i][j + 1])->getOutVEGF() +
                    ((OxyCell *)m_map[0][i + 1][j])->getOutVEGF() -
                    4.0 * ((OxyCell *)m_map[0][i][j])->getOutVEGF());
        }
    }

    for(int l(0); l < m_nlayer; l++){
        for(int i(0); i < m_nrow; i++){
            for(int j(0); j < m_ncol; j++){
                ((OxyCell *)m_map[l][i][j])->setInDiffO2(diffO2[l][i][j]);
                ((OxyCell *)m_map[l][i][j])->setInDiffVEGF(diffVegf[l][i][j]);
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
                if(fabs(((OxyCell *)m_map[l][i][j])->getOutPO2() -
                        ((OxyCell *)m_map[l][i][j])->getPO2()) >
                        1e-2 * ((OxyCell *)m_map[l][i][j])->getOutPO2() ||
                        fabs(((OxyCell *)m_map[l][i][j])->getOutVEGF() -
                             ((OxyCell *)m_map[l][i][j])->getVEGF()) >
                        1e-2 * ((OxyCell *)m_map[l][i][j])->getOutVEGF()){
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
