/**
 * @file diffusion3D.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef DEF_OXYTISSUE
#define DEF_OXYTISSUE

#include <string>
#include <vector>

#include "model.hpp"
#include "oxyCell.hpp"

//State variables
#define ST_OXYSTABLE          m_st->at(0)
#define ST_VEGFSTABLE         m_st->at(1)
#define ST_TIME_TO_OXYSTABLE  m_st->at(2)
#define ST_TIME_TO_VEGFSTABLE m_st->at(3)

//Outputs
#define OUT_HYP_DENS           m_out->at(0)
#define OUT_PO2_MED            m_out->at(1)
#define OUT_PO2_MEAN           m_out->at(2)
#define OUT_VEGF_MED           m_out->at(3)
#define OUT_VEGF_MEAN          m_out->at(4)
#define OUT_OXYSTABLE          m_out->at(5)
#define OUT_VEGFSTABLE         m_out->at(6)
#define OUT_TIME_TO_OXYSTABLE  m_out->at(7)
#define OUT_TIME_TO_VEGFSTABLE m_out->at(8)

//Internal parameters
#define PAR_OXY_ANG m_param->at(0)
#define PAR_DVEGF   m_param->at(1)
#define PAR_DO2	    m_param->at(2)

class OxyTissue : public Model{
public:
    OxyTissue(const int nrow, const int ncol, const int nlayer,
              const double cellSize, const std::string nFInVes,
              const double ang, const double Dvegf, const double VmaxVegf,
              const double KmVegf, const double hypVegf, const double DO2,
              const double VmaxO2, const double KmO2, const double pO2NormVes,
              const double pO2TumVes, const double hypThres) ;
    OxyTissue(const int nrow, const int ncol, const int nlayer,
              const double cellSize, const std::vector<bool> &inVes,
              const double ang, const double Dvegf, const double VmaxVegf,
              const double KmVegf, const double hypVegf, const double DO2,
              const double VmaxO2, const double KmO2, const double pO2NormVes,
              const double pO2TumVes, const double hypThres);
    virtual ~OxyTissue();
    virtual int initModel();
    virtual int calcModelOut();
    virtual int updateModel(double currentTime, const double DT);
    virtual int terminateModel();
    int getNumHyp() const;
    bool isOxyStable() const;
    bool isVegfStable() const;

protected:
    int m_ncol, m_nlayer, m_nrow;
    OxyCell ****m_map;
};

#endif
