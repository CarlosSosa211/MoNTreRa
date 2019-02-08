/**
 * @file prostateCell.hpp
 * @brief
 * @author Nicolas Ciferri
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef DEF_OXYCELL
#define DEF_OXYCELL

#include <cmath>

#include "absOxyCell.hpp"

//Inputs
#define IN_DIFF_O2   m_in->at(3)
#define IN_CONS_O2   m_in->at(4)
#define IN_DIFF_VEGF m_in->at(5)
#define IN_CONS_VEGF m_in->at(6)

//Internal parameters
#define PAR_OXYCELL_ANG m_param->at(5)
#define PAR_VMAX_VEGF   m_param->at(6)
#define PAR_KM_VEGF     m_param->at(7)
#define PAR_VMAX_O2     m_param->at(8)
#define PAR_KM_O2       m_param->at(9)



class OxyCell: public AbsOxyCell{
public :
    OxyCell();
    OxyCell(const double ang, const double  VmaxVegf, const double KmVegf,
             const double hypVegf, const double VmaxO2, const double KmO2,
             const double pO2NormVes, const double pO2TumVes,
             const double hypThres, Model *const parent);
    virtual ~OxyCell();
    virtual int updateModel(const double currentTime, const double DT);
    void addToEdge(OxyCell *const cell);
    void calcConsO2();
    void calcConsVegf();
    std::vector<OxyCell *> *getEdge() const;
    void setInDiffO2(const double input);
    void setInDiffVEGF(const double input);

protected:
    std::vector<OxyCell *> *m_edge;
};

#endif
