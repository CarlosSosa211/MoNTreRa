/**
 * @file absOxyCell.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.07.17
 */

#ifndef DEF_ABSTOXYCELL
#define DEF_ABSTOXYCELL

#include "model.hpp"

//Inputs
#define IN_OXYDEAD     m_in->at(0)
#define IN_OXYNORM_VES m_in->at(1)
#define IN_OXYTUM_VES  m_in->at(2)

//State variables
#define ST_PO2         m_st->at(0)
#define ST_VEGF        m_st->at(1)
#define ST_HYP         m_st->at(2)
#define ST_OXYDEAD     m_st->at(3)
#define ST_OXYNORM_VES m_st->at(4)
#define ST_OXYTUM_VES  m_st->at(5)

//Outputs
#define OUT_PO2       m_out->at(0)
#define OUT_VEGF      m_out->at(1)

//Internal parameters
#define PAR_OXYPO2       m_param->at(0)
#define PAR_PO2_NORM_VES m_param->at(1)
#define PAR_PO2_TUM_VES  m_param->at(2)
#define PAR_HYP_THRES    m_param->at(3)
#define PAR_HYP_VEGF     m_param->at(4)

class AbsOxyCell: public Model{
public :
    AbsOxyCell();
    virtual ~AbsOxyCell();
    virtual int calcModelOut();
    virtual int initModel();
    virtual int terminateModel();
    virtual int updateModel(const double currentTime, const double DT) = 0;
    bool getHyp() const;
    double getOutPO2() const;
    double getOutVEGF() const;
    double getPO2() const;
    double getVEGF() const;
    void setInDead(const double input);
    void setInNormVes(const double input);
    void setInTumVes(const double input);
};

#endif
