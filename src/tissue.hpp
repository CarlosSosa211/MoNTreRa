/**
 * @file tissue.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef DEF_TISSUE
#define DEF_TISSUE

#include <string>
#include <vector>

#include "cell.hpp"
#include "model.hpp"
#include "treatment.hpp"

//State variables
#define ST_TUM_DENS       m_st->at(0)
#define ST_PREV_TUM_DENS  m_st->at(1)
#define ST_INT_TUM_DENS   m_st->at(2)

//Outputs
#define OUT_TUM_DENS      m_out->at(0)
#define OUT_G1_DENS       m_out->at(1)
#define OUT_S_DENS        m_out->at(2)
#define OUT_G2_DENS       m_out->at(3)
#define OUT_M_DENS        m_out->at(4)
#define OUT_G0_DENS       m_out->at(5)
#define OUT_VES_DENS      m_out->at(6)
#define OUT_NORM_VES_DENS m_out->at(7)
#define OUT_TUM_VES_DENS  m_out->at(8)
#define OUT_TIME_TO_50    m_out->at(9)
#define OUT_TIME_TO_80    m_out->at(10)
#define OUT_TIME_TO_90    m_out->at(11)
#define OUT_TIME_TO_95    m_out->at(12)
#define OUT_TIME_TO_99    m_out->at(13)
#define OUT_TIME_TO_999   m_out->at(14)
#define OUT_DOSE_TO_50    m_out->at(15)
#define OUT_DOSE_TO_80    m_out->at(16)
#define OUT_DOSE_TO_90    m_out->at(17)
#define OUT_DOSE_TO_95    m_out->at(18)
#define OUT_DOSE_TO_99    m_out->at(19)
#define OUT_DOSE_TO_999   m_out->at(20)
#define OUT_KILLED_CELLS  m_out->at(21)
#define OUT_INT_TUM_DENS  m_out->at(22)

//Internal parameters
#define PAR_INIT_TUM_DENS m_param->at(0)
#define PAR_INIT_VES_DENS m_param->at(1)


class Tissue : public Model{
public:
    Tissue(const int nrow, const int ncol,
           const int nlayer, Treatment *const treatment);
    Tissue(const int nrow, const int ncol, const int nlayer,
           const std::string nFInTum, const std::string nFInVes,
           const double tumGrowth, const double doubTime, const int edgeOrder,
           std::vector<double> cycDur, std::vector<double> cycDistrib,
           const double res, const double fibDoubTime, const double ang,
           const double angTime, const double vegfThres,
           std::vector<double> alpha, std::vector<double> beta,
           const double doseThres, const double arresTime, 
           Treatment *const treatment, const double hypNecThres);
    Tissue(const int nrow, const int ncol, const int nlayer,
           const std::vector<bool> &inTum, const std::vector<bool> &inVes,
           const double tumGrowth, const double doubTime, const int edgeOrder,
           std::vector<double> cycDur, std::vector<double> cycDistrib,
           const double res, const double fibDoubTime, const double ang,
           const double angTime, const double vegfThres,
           std::vector<double> alpha, std::vector<double> beta,
           const double doseThres, const double arrestTime, 
           Treatment *const treatment, const double hypNecThres);
    virtual ~Tissue();
    virtual int calcModelOut();
    virtual int initModel();
    virtual int terminateModel();
    virtual int updateModel(const double currentTime,
                            const double DT);
    int getNumFib() const;
    int getNumDead() const;
    int getNumG0() const;
    int getNumG1() const;
    int getNumG2() const;
    int getNumM() const;
    int getNumNormVes() const;
    int getNumS() const;
    int getNumTum() const;
    int getNumTumVes() const;
    int getNumVes() const;
    Treatment *getTreatment() const;
    void printNeededDoseAndTime(int sel) const;

protected:
    int m_ncol, m_nlayer, m_nrow;
    double m_doseNeeded[6], m_timeNeeded[6];
    Treatment *m_treatment;
};

#endif
