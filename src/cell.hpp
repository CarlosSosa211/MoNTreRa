/**
 * @file cell.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef DEF_CELL
#define DEF_CELL

#include <vector>

#include "model.hpp"
#include "tissue.hpp"
#include "treatment.hpp"

//Inputs
#define IN_FIB       m_in->at(0)
#define IN_TUM       m_in->at(1)
#define IN_NORM_VES  m_in->at(2)
#define IN_TUM_VES   m_in->at(3)
#define IN_HYP_NEC   m_in->at(4)
#define IN_MIT_CAT   m_in->at(5)
#define IN_APOP      m_in->at(6)
#define IN_TIMER     m_in->at(7)
#define IN_PO2       m_in->at(8)
#define IN_VEGF      m_in->at(9)

//State variables
#define ST_FIB      m_st->at(0)
#define ST_TUM      m_st->at(1)
#define ST_NORM_VES m_st->at(2)
#define ST_TUM_VES  m_st->at(3)
#define ST_VES      m_st->at(4)
#define ST_HYP_NEC  m_st->at(5)
#define ST_MIT_CAT  m_st->at(6)
#define ST_APOP     m_st->at(7)
#define ST_DEAD     m_st->at(8)
#define ST_G1       m_st->at(9)
#define ST_S        m_st->at(10)
#define ST_G2       m_st->at(11)
#define ST_M        m_st->at(12)
#define ST_G0       m_st->at(13)
#define ST_TIMER    m_st->at(14)
#define ST_VEGF     m_st->at(15)
#define ST_DAM      m_st->at(16)
#define ST_ARREST   m_st->at(17)
#define ST_ACC_DOSE m_st->at(18)
#define ST_PO2      m_st->at(19)

//Outputs
#define OUT_STATE m_out->at(0)
#define OUT_CYCLE m_out->at(1)

//Internal parameters
#define PAR_TUM_GROWTH     m_param->at(0)
#define PAR_DOUB_TIME      m_param->at(1)
#define PAR_LIM_G1S        m_param->at(2)
#define PAR_LIM_SG2        m_param->at(3)
#define PAR_LIM_G2M        m_param->at(4)
#define PAR_RES            m_param->at(5)
#define PAR_FIB_DOUB_TIME  m_param->at(6)
#define PAR_ANG            m_param->at(7)
#define PAR_ANG_TIME       m_param->at(8)
#define PAR_VEGF_THRES     m_param->at(9)
#define PAR_ALPHA          m_param->at(10)
#define PAR_ALPHA_FIB      m_param->at(11)
#define PAR_ALPHA_G1       m_param->at(12)
#define PAR_ALPHA_S        m_param->at(13)
#define PAR_ALPHA_G2       m_param->at(14)
#define PAR_ALPHA_M        m_param->at(15)
#define PAR_ALPHA_G0       m_param->at(16)
#define PAR_ALPHA_NORM_VES m_param->at(17)
#define PAR_ALPHA_TUM_VES  m_param->at(18)
#define PAR_BETA           m_param->at(19)
#define PAR_BETA_FIB       m_param->at(20)
#define PAR_BETA_G1        m_param->at(21)
#define PAR_BETA_S         m_param->at(22)
#define PAR_BETA_G2        m_param->at(23)
#define PAR_BETA_M         m_param->at(24)
#define PAR_BETA_G0        m_param->at(25)
#define PAR_BETA_NORM_VES  m_param->at(26)
#define PAR_BETA_TUM_VES   m_param->at(27)
#define PAR_DOSE_THRES     m_param->at(28)
#define PAR_ARREST_TIME    m_param->at(29)
#define PAR_M              m_param->at(30)
#define PAR_K              m_param->at(31)
#define PAR_OXY            m_param->at(32)
#define PAR_HYP_NEC_THRES  m_param->at(33)

class Cell: public Model{
public :
    Cell(Model *const parent = 0);
    Cell(const int i, const int j, const int l,
         const double tumGrowth, const double doubTime,
         std::vector<double> cycDur, const double res,
         const double fibDoubTime, const double ang,
         const double angTime, const double vegfThres,
         std::vector<double> alpha, std::vector<double> beta,
         const double doseThres, const double arrestTime, 
         const double oxy, const double hypNecThres,
         Model *const parent = 0);
    virtual ~Cell();
    virtual int calcModelOut();
    virtual int initModel();
    virtual int terminateModel();
    virtual int updateModel(const double currentTime,
                            const double DT);
    void addToEdge(Cell *const cell);
    void calcFibProlif(double DT);
    double calcG1SF() const;
    double calcG2MF() const;
    void calcHypNec();
    void calcNormVesProlif(double DT);
    double calcOER() const;
    void calcRespToIrr();
    double calcSF() const;
    void calcTumGrowth(double DT);
    void calcTumVesProlif(double DT);
    bool getFib() const;
    bool getApop() const;
    double getAccDose() const;
    bool getDead() const;
    double getDoubTime() const;
    std::vector<Cell *> *getEdge() const;
    bool getG0() const;
    bool getG1() const;
    bool getG2() const;
    bool getHypNec() const;
    bool getM() const;
    bool getMitCat() const;
    bool getNormVes() const;
    int getOutState() const;
    double getR() const;
    bool getS() const;
    bool getTum() const;
    bool getTumDam() const;
    bool getTumNotDam() const;
    bool getTumVes() const;
    bool getVes() const;
    Cell *searchInitSpaceForTum() const;
    Cell *searchSpaceForFib() const;
    Cell *searchSpaceForTum() const;
    Cell *searchSpaceForTumVes() const;
    Cell *searchSpaceForVes() const;
    void setInFib(const double input);
    void setInApop(const double input);
    void setInHypNec(const double input);
    void setInMitCat(const double input);
    void setInNormVes(const double input);
    void setInPO2(const double input);
    void setInTimer(const double input);
    void setInTum(const double input);
    void setInTumVes(const double input);
    void setInVegf(const double input);
    void setR(const Cell *origCell);

protected:
    int m_i, m_j, m_l, m_r;
    Treatment *m_treatment;
    std::vector<Cell *> *m_edge;
};

bool compRedge(Cell *a, Cell *b);
#endif
