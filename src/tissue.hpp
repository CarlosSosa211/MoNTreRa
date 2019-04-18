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
#define ST_TUM_DENS           m_st->at(0)
#define ST_PREV_TUM_DENS      m_st->at(1)
#define ST_INT_TUM_DENS       m_st->at(2)
#define ST_END_TREAT_TUM_DENS m_st->at(3)
#define ST_3MON_TUM_DENS      m_st->at(4)
#define ST_REC                m_st->at(5)
#define ST_COUNT_REC          m_st->at(6)
#define ST_REC_TUM_DENS       m_st->at(7)
#define ST_REC_TIME           m_st->at(8)
#define ST_50_KILLED          m_st->at(9)
#define ST_80_KILLED          m_st->at(10)
#define ST_90_KILLED          m_st->at(11)
#define ST_95_KILLED          m_st->at(12)
#define ST_99_KILLED          m_st->at(13)
#define ST_999_KILLED         m_st->at(14)
#define ST_TIME_TO_50         m_st->at(15)
#define ST_TIME_TO_80         m_st->at(16)
#define ST_TIME_TO_90         m_st->at(17)
#define ST_TIME_TO_95         m_st->at(18)
#define ST_TIME_TO_99         m_st->at(19)
#define ST_TIME_TO_999        m_st->at(20)
#define ST_DOSE_TO_50         m_st->at(21)
#define ST_DOSE_TO_80         m_st->at(22)
#define ST_DOSE_TO_90         m_st->at(23)
#define ST_DOSE_TO_95         m_st->at(24)
#define ST_DOSE_TO_99         m_st->at(25)
#define ST_DOSE_TO_999        m_st->at(26)
#define ST_CONTROLLED         m_st->at(27)
#define ST_DOSE_TO_CONTROL    m_st->at(28)
#define ST_VES_DENS           m_st->at(29)
#define ST_NORM_VES_DENS      m_st->at(30)
#define ST_TUM_VES_DENS       m_st->at(31)
#define ST_DEAD_DENS          m_st->at(32)

//Outputs
#define OUT_TUM_DENS           m_out->at(0)
#define OUT_G1_DENS            m_out->at(1)
#define OUT_S_DENS             m_out->at(2)
#define OUT_G2_DENS            m_out->at(3)
#define OUT_M_DENS             m_out->at(4)
#define OUT_G0_DENS            m_out->at(5)
#define OUT_VES_DENS           m_out->at(6)
#define OUT_NORM_VES_DENS      m_out->at(7)
#define OUT_TUM_VES_DENS       m_out->at(8)
#define OUT_TIME_TO_50         m_out->at(9)
#define OUT_TIME_TO_80         m_out->at(10)
#define OUT_TIME_TO_90         m_out->at(11)
#define OUT_TIME_TO_95         m_out->at(12)
#define OUT_TIME_TO_99         m_out->at(13)
#define OUT_TIME_TO_999        m_out->at(14)
#define OUT_DOSE_TO_50         m_out->at(15)
#define OUT_DOSE_TO_80         m_out->at(16)
#define OUT_DOSE_TO_90         m_out->at(17)
#define OUT_DOSE_TO_95         m_out->at(18)
#define OUT_DOSE_TO_99         m_out->at(19)
#define OUT_DOSE_TO_999        m_out->at(20)
#define OUT_KILLED_CELLS       m_out->at(21)
#define OUT_INT_TUM_DENS       m_out->at(22)
#define OUT_TUM_VOL            m_out->at(23)
#define OUT_END_TREAT_TUM_DENS m_out->at(24)
#define OUT_3MON_TUM_DENS      m_out->at(25)
#define OUT_REC_TUM_DENS       m_out->at(26)
#define OUT_REC_TIME           m_out->at(27)
#define OUT_50_KILLED          m_out->at(28)
#define OUT_80_KILLED          m_out->at(29)
#define OUT_90_KILLED          m_out->at(30)
#define OUT_95_KILLED          m_out->at(31)
#define OUT_99_KILLED          m_out->at(32)
#define OUT_999_KILLED         m_out->at(33)
#define OUT_REC                m_out->at(34)
#define OUT_CONTROLLED         m_out->at(35)
#define OUT_DOSE_TO_CONTROL    m_out->at(36)
#define OUT_DEAD_DENS          m_out->at(37)

//Internal parameters
#define PAR_INIT_TUM_DENS m_param->at(0)
#define PAR_INIT_VES_DENS m_param->at(1)

class Tissue : public Model{
public:
    Tissue(const int nrow, const int ncol, const int nlayer,
           Treatment *const treatment);
    Tissue(const int nrow, const int ncol, const int nlayer,
           const double cellSize, const std::string nFInTum,
           const std::string nFInVes, const bool tumGrowth,
           const double doubTime, const int edgeOrder,
           std::vector<double> cycDur, std::vector<double> cycDistrib,
           const bool res, const double fibDoubTime, const bool ang,
           const double angTime, const double vegfThres,
           std::vector<double> alpha, std::vector<double> beta,
           Treatment *const treatment, const double doseThres,
           const double arresTime, const int oxy, const double hypNecThres);
    Tissue(const int nrow, const int ncol, const int nlayer,
           const double cellSize, const std::vector<bool> &inTum,
           const std::vector<bool> &inVes, const bool tumGrowth,
           const double doubTime, const int edgeOrder,
           std::vector<double> cycDur, std::vector<double> cycDistrib,
           const bool res, const double fibDoubTime, const bool ang,
           const double angTime, const double vegfThres,
           std::vector<double> alpha, std::vector<double> beta,
           Treatment *const treatment, const double doseThres,
           const double arrestTime, const int oxy, const double hypNecThres);
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
    int getNumTumNotDam() const;
    int getNumTumVes() const;
    int getNumVes() const;
    Treatment *getTreatment() const;

protected:
    int m_ncol, m_nlayer, m_nrow;
    double m_cellSize;
    Treatment *m_treatment;
};

#endif
