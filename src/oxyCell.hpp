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

#include "absOxyCell.hpp"

//Inputs
#define IN_DIFF_O2   m_in->at(3)
#define IN_CONS_O2   m_in->at(4)
#define IN_DIFF_VEGF m_in->at(5)
#define IN_CONS_VEGF m_in->at(6)

//Internal parameters
#define PAR_VMAX      m_param->at(5)
#define PAR_KM        m_param->at(6)
#define PAR_VMAX_VEGF m_param->at(7)
#define PAR_KM_VEGF   m_param->at(8)


class OxyCell: public AbsOxyCell{
public :
    OxyCell();
    OxyCell(const double Vmax, const double Km,
            const double pO2NormVes, const double pO2TumVes,
            const double hypThres, const double  VmaxVegf,
            const double KmVegf, const double hypVegf,
            Model *const parent);
    virtual ~OxyCell();
    virtual int updateModel(const double currentTime, const double DT);
    void calcConsO2();
    void calcConsVegf();
    void setInDiffO2(const double input);
    void setInDiffVEGF(const double input);
};

#endif
