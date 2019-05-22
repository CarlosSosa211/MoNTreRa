/**
 * @file diffusion3D.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.19.17
 */

#ifndef DEF_OXYTISSUE
#define DEF_OXYTISSUE

#include <string>
#include <vector>

#include "absOxyTissue.hpp"
#include "model.hpp"
#include "oxyCell.hpp"

//Inputs
#define OXYTISSUE_NUM_IN_B 0
#define OXYTISSUE_NUM_IN_I 0
#define OXYTISSUE_NUM_IN_D 0

//State variables
#define OXYTISSUE_NUM_ST_B 2
#define OXYTISSUE_NUM_ST_I 0
#define OXYTISSUE_NUM_ST_D 2

#define ST_OXY_STABLE  m_stB[0] //pO2 stability
#define ST_VEGF_STABLE m_stB[1] //VEGF concentration stability

#define ST_TIME_TO_OXY_STABLE  m_stD[0] //time to pO2 stability (ms)
#define ST_TIME_TO_VEGF_STABLE m_stD[1] //time to VEGF conc. stability (ms)

//Outputs
#define OXYTISSUE_NUM_OUT_B 2
#define OXYTISSUE_NUM_OUT_I 0
#define OXYTISSUE_NUM_OUT_D 7

#define OUT_OXYSTABLE  m_outB[0] //pO2 stability
#define OUT_VEGFSTABLE m_outB[1] //VEGF stability

#define OUT_TIME_TO_OXYSTABLE  m_outD[5] //time to pO2 stability (ms)
#define OUT_TIME_TO_VEGFSTABLE m_outD[6] //time to VEGF conc. stability (ms)

//Internal parameters
#define OXYTISSUE_NUM_PAR_B 1
#define OXYTISSUE_NUM_PAR_I 1
#define OXYTISSUE_NUM_PAR_D 2

#define PAR_DVEGF m_parD[0] //VEGF diffusion coefficient (um^2/ms)
#define PAR_DO2	  m_parD[1] //pO2 diffusion coefficient (um^2/ms)


class OxyTissue : public AbsOxyTissue{
public:
    OxyTissue(const int nrow, const int ncol, const int nlayer,
              const double cellSize, const std::vector<bool> &inVes,
              const bool ang, const double DVegf, const double VmaxVegf,
              const double KmVegf, const double hypVegf, const int oxy,
              const double DO2, const double VmaxO2, const double KmO2,
              const double pO2NormVes, const double pO2TumVes,
              const double hypThres);
    virtual ~OxyTissue();
    virtual int initModel();
    virtual int calcModelOut();
    virtual int updateModel(double currentTime, const double DT);
    bool isOxyStable() const;
    bool isVegfStable() const;

protected:
    int m_ncol, m_nlayer, m_nrow;
    OxyCell ****m_map;
};

#endif
