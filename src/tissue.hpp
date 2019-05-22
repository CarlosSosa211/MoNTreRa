/**
 * @file tissue.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.19.17
 */

#ifndef DEF_TISSUE
#define DEF_TISSUE

#include <string>
#include <vector>

#include "cell.hpp"
#include "model.hpp"
#include "treatment.hpp"

//Inputs
#define TISSUE_NUM_IN_B 0
#define TISSUE_NUM_IN_I 0
#define TISSUE_NUM_IN_D 0

//State variables
#define TISSUE_NUM_ST_B 8
#define TISSUE_NUM_ST_I 1
#define TISSUE_NUM_ST_D 24

#define ST_REC        m_stB[0] //recurrence
#define ST_50_KILLED  m_stB[1] //50% of initial tumour cells killed
#define ST_80_KILLED  m_stB[2] //80% of initial tumour cells killed
#define ST_90_KILLED  m_stB[3] //90% of initial tumour cells killed
#define ST_95_KILLED  m_stB[4] //95% of initial tumour cells killed
#define ST_99_KILLED  m_stB[5] //99% of initial tumour cells killed
#define ST_999_KILLED m_stB[6] //99.9% of initial tumour cells killed
#define ST_CONTROLLED m_stB[7] //tumour controlled

#define ST_COUNT_REC m_stI[0] //count of tum. growth iter. to rec.

#define ST_TUM_DENS           m_stD[0]  //tumour density
#define ST_PREV_TUM_DENS      m_stD[1]  //previous tumour density
#define ST_INT_TUM_DENS       m_stD[2]  //integral of tumour density
#define ST_END_TREAT_TUM_DENS m_stD[3]  //tumour density at the end of treatment
#define ST_3MON_TUM_DENS      m_stD[4]  //tum. dens. 3 M after start of treat.
#define ST_REC_TUM_DENS       m_stD[5]  //tum. dens. at start of recurrence
#define ST_REC_TIME           m_stD[6]  //time at start of recurrence (h)
#define ST_TIME_TO_50         m_stD[7]  //t to kill 50% of init. tum. cells (h)
#define ST_TIME_TO_80         m_stD[8]  //t to kill 50% of init. tum. cells (h)
#define ST_TIME_TO_90         m_stD[9] //t to kill 90% of init. tum. cells (h)
#define ST_TIME_TO_95         m_stD[10] //t to kill 95% of init. tum. cells (h)
#define ST_TIME_TO_99         m_stD[11] //t to kill 99% of init. tum. cells (h)
#define ST_TIME_TO_999        m_stD[12] //t to kill 99.9% of in. tum. cells (h)
#define ST_DOSE_TO_50         m_stD[13] //d to kill 50% of init. tum. cells (Gy)
#define ST_DOSE_TO_80         m_stD[14] //d to kill 80% of init. tum. cells (Gy)
#define ST_DOSE_TO_90         m_stD[15] //d to kill 90% of init. tum. cells (Gy)
#define ST_DOSE_TO_95         m_stD[16] //d to kill 95% of init. tum. cells (Gy)
#define ST_DOSE_TO_99         m_stD[17] //d to kill 99% of init. tum. cells (Gy)
#define ST_DOSE_TO_999        m_stD[18] //d to kill 99.9% of in. tum. cells (Gy)
#define ST_DOSE_TO_CONTROL    m_stD[19] //d to tumour control (Gy)
#define ST_VES_DENS           m_stD[20] //vascular density
#define ST_NORM_VES_DENS      m_stD[21] //pre-existing vascular density
#define ST_TUM_VES_DENS       m_stD[22] //neo-created vascular density
#define ST_DEAD_DENS          m_stD[23] //dead cell density

//Outputs
#define TISSUE_NUM_OUT_B 8
#define TISSUE_NUM_OUT_I 0
#define TISSUE_NUM_OUT_D 30

#define OUT_50_KILLED  m_outB[0] //50% of initial tumour cells killed
#define OUT_80_KILLED  m_outB[1] //80% of initial tumour cells killed
#define OUT_90_KILLED  m_outB[2] //90% of initial tumour cells killed
#define OUT_95_KILLED  m_outB[3] //95% of initial tumour cells killed
#define OUT_99_KILLED  m_outB[4] //99% of initial tumour cells killed
#define OUT_999_KILLED m_outB[5] //99.9% of initial tumour cells killed
#define OUT_REC        m_outB[6] //recurrence
#define OUT_CONTROLLED m_outB[7] //tumour controlled

#define OUT_END_TREAT_TUM_DENS m_outD[0]  //tum. dens. at the end of treatment
#define OUT_3MON_TUM_DENS      m_outD[1]  //tum. dens. 3 M after start of treat.
#define OUT_INT_TUM_DENS       m_outD[2]  //integral of tumour density
#define OUT_TIME_TO_50         m_outD[3]  //t to kill 50% of init. tum. cells
#define OUT_TIME_TO_80         m_outD[4]  //t to kill 80% of init. tum. cells
#define OUT_TIME_TO_90         m_outD[5]  //t to kill 90% of init. tum. cells
#define OUT_TIME_TO_95         m_outD[6]  //t to kill 95% of init. tum. cells
#define OUT_TIME_TO_99         m_outD[7]  //t to kill 99% of init. tum. cells
#define OUT_TIME_TO_999        m_outD[8]  //t to kill 99.9% of init. tum. cells
#define OUT_DOSE_TO_50         m_outD[9]  //d to kill 50% of init. tum. cells
#define OUT_DOSE_TO_80         m_outD[10] //d to kill 80% of init. tum. cells
#define OUT_DOSE_TO_90         m_outD[11] //d to kill 90% of init. tum. cells
#define OUT_DOSE_TO_95         m_outD[12] //d to kill 95% of init. tum. cells
#define OUT_DOSE_TO_99         m_outD[13] //d to kill 99% of init. tum. cells
#define OUT_DOSE_TO_999        m_outD[14] //d to kill 99.9% of init. tum. cells
#define OUT_DOSE_TO_CONTROL    m_outD[15] //dose to tumour control (Gy)
#define OUT_REC_TUM_DENS       m_outD[16] //tum. dens. at start of recurrence
#define OUT_REC_TIME           m_outD[17] //time at start of recurrence (h)
#define OUT_TUM_DENS           m_outD[18] //tumour density
#define OUT_TUM_VOL            m_outD[19] //tumour volume (mm^3)
#define OUT_VES_DENS           m_outD[20] //vascular density
#define OUT_NORM_VES_DENS      m_outD[21] //pre-existing vascular density
#define OUT_TUM_VES_DENS       m_outD[22] //neo-created vascular density
#define OUT_KILLED_CELLS       m_outD[23] //killed tumour cells percentage
#define OUT_DEAD_DENS          m_outD[24] //dead cell density
#define OUT_G1_DENS            m_outD[25] //G1 tumour cells percentage
#define OUT_S_DENS             m_outD[26] //S tumour cells percentage
#define OUT_G2_DENS            m_outD[27] //G2 tumour cells percentage
#define OUT_M_DENS             m_outD[28] //M tumour cells percentage
#define OUT_G0_DENS            m_outD[29] //G0 tumour cells percentage

//Internal parameters
#define TISSUE_NUM_PAR_B 0
#define TISSUE_NUM_PAR_I 0
#define TISSUE_NUM_PAR_D 2

#define PAR_INIT_TUM_DENS m_parD[0] //intial tumour density
#define PAR_INIT_VES_DENS m_parD[1] //initial vascular density

class Tissue : public Model{
public:
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
    virtual int updateModel(const double currentTime, const double DT);
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
