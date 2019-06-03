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
#define OXYCELL_NUM_IN_B 3
#define OXYCELL_NUM_IN_I 0
#define OXYCELL_NUM_IN_D 4

#define IN_DIFF_O2   m_inD[0] //diffused O2
#define IN_CONS_O2   m_inD[1] //consumed O2
#define IN_DIFF_VEGF m_inD[2] //diffused VEGF
#define IN_CONS_VEGF m_inD[3] //consumed O2

//State variables
#define OXYCELL_NUM_ST_B 6
#define OXYCELL_NUM_ST_I 0
#define OXYCELL_NUM_ST_D 2

//Outputs
#define OXYCELL_NUM_OUT_B 0
#define OXYCELL_NUM_OUT_I 0
#define OXYCELL_NUM_OUT_D 2

//Internal parameters
#define OXYCELL_NUM_PAR_B 1
#define OXYCELL_NUM_PAR_I 0
#define OXYCELL_NUM_PAR_D 8

#define PAR_OXYCELL_ANG m_parB[0] //angiogenesis

#define PAR_VMAX_VEGF   m_parD[4] //maximum VEGF consumption ratio (mol/um^3ms)
#define PAR_KM_VEGF     m_parD[5] //VEGF Michaelis constant (mol/um^3)
#define PAR_VMAX_O2     m_parD[6] //maximum pO2 consumption ratio (mmHg/ms)
#define PAR_KM_O2       m_parD[7] //pO2 Michaelis constant (mmHg)


class OxyCell: public AbsOxyCell{
public :
    OxyCell();
    OxyCell(const bool ang, const double  VmaxVegf, const double KmVegf,
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
