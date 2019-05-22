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
#define CELL_NUM_IN_B 7
#define CELL_NUM_IN_I 0
#define CELL_NUM_IN_D 3

#define IN_FIB      m_inB[0] //healthy cell
#define IN_TUM      m_inB[1] //tumour cell
#define IN_NORM_VES m_inB[2] //pre-existing endothelial cell
#define IN_TUM_VES  m_inB[3] //neo-created endothelial cell
#define IN_HYP_NEC  m_inB[4] //hypoxic necrotic cell
#define IN_MIT_CAT  m_inB[5] //mitotic dead cell
#define IN_APOP     m_inB[6] //apoptotic cell

#define IN_TIMER m_inD[0]
#define IN_PO2   m_inD[1]
#define IN_VEGF  m_inD[2]

//State variables
#define CELL_NUM_ST_B 15
#define CELL_NUM_ST_I 0
#define CELL_NUM_ST_D 7

#define ST_FIB      m_stB[0]  //healthy cell
#define ST_TUM      m_stB[1]  //tumour cell
#define ST_NORM_VES m_stB[2]  //pre-existing endothelial cell
#define ST_TUM_VES  m_stB[3]  //neo-created endothelial cell
#define ST_VES      m_stB[4]  //endothelial cell
#define ST_HYP_NEC  m_stB[5]  //hypoxic necrotic cell
#define ST_MIT_CAT  m_stB[6]  //mitotic dead cell
#define ST_APOP     m_stB[7]  //apoptotic cell
#define ST_DEAD     m_stB[8]  //dead cell
#define ST_G1       m_stB[9]  //G1 phase cell
#define ST_S        m_stB[10] //S phase cell
#define ST_G2       m_stB[11] //G2 phase cell
#define ST_M        m_stB[12] //M phase cell
#define ST_G0       m_stB[13] //G0 phase cell
#define ST_DAM      m_stB[14] //damaged cell

#define ST_TIMER    m_stD[0] //timer (h)
#define ST_VEGF     m_stD[1] //VEGF concentration
#define ST_ALPHA    m_stD[2] //alpha (Gy^-1)
#define ST_BETA     m_stD[3] //beta (Gy^-2)
#define ST_ARREST   m_stD[4] //arrest time (h)
#define ST_ACC_DOSE m_stD[5] //accumulated dose (Gy)
#define ST_PO2      m_stD[6] //pO2 (mmHg)

//Outputs
#define CELL_NUM_OUT_B 15
#define CELL_NUM_OUT_I 0
#define CELL_NUM_OUT_D 5

#define OUT_STATE m_outI[0] //type of cell
#define OUT_CYCLE m_outI[1] //phase of the cell

//Internal parameters
#define CELL_NUM_PAR_B 3
#define CELL_NUM_PAR_I 1
#define CELL_NUM_PAR_D 28

#define PAR_TUM_GROWTH m_parB[0] //tumour growth
#define PAR_RES        m_parB[1] //healthy cell division
#define PAR_ANG        m_parB[2] //angiogenesis

#define PAR_OXY m_parI[0]

#define PAR_DOUB_TIME      m_parD[0]  //dur. of the cycle of tumour cells (h)
#define PAR_LIM_G1S        m_parD[1]  //G1/S transition time (h)
#define PAR_LIM_SG2        m_parD[2]  //S/G2 transition time (h)
#define PAR_LIM_G2M        m_parD[3]  //G2/M transition time (h)
#define PAR_FIB_DOUB_TIME  m_parD[4]  //dur. of the cycle of healthy cells (h)
#define PAR_ANG_TIME       m_parD[5]  //dur. of the cycle of end. cells (h)
#define PAR_VEGF_THRES     m_parD[6]  //VEGF thres. to end. cell div. (mol/um^3)
#define PAR_ALPHA_FIB      m_parD[7]  //alpha of healthy cells (Gy^-1)
#define PAR_ALPHA_G1       m_parD[8]  //alpha of G1 tumour cells (Gy^-1)
#define PAR_ALPHA_S        m_parD[9]  //alpha of S tumour cells (Gy^-1)
#define PAR_ALPHA_G2       m_parD[10] //alpha of G2 tumour cells (Gy^-1)
#define PAR_ALPHA_M        m_parD[11] //alpha of M tumour cells (Gy^-1)
#define PAR_ALPHA_G0       m_parD[12] //alpha of G0 tumour cells (Gy^-1)
#define PAR_ALPHA_NORM_VES m_parD[13] //alpha of pre-exist. end. cells (Gy^-1)
#define PAR_ALPHA_TUM_VES  m_parD[14] //alpha of neo-crea. end. cells (Gy^-1)
#define PAR_BETA_FIB       m_parD[15] //beta of healthy cells (Gy^-2)
#define PAR_BETA_G1        m_parD[16] //beta of G1 tumour cells (Gy^-2)
#define PAR_BETA_S         m_parD[17] //beta of S tumour cells (Gy^-2)
#define PAR_BETA_G2        m_parD[18] //beta of G2 tumour cells (Gy^-2)
#define PAR_BETA_M         m_parD[19] //beta of M tumour cells (Gy^-2)
#define PAR_BETA_G0        m_parD[20] //beta of G0 tumour cells (Gy^-2)
#define PAR_BETA_NORM_VES  m_parD[21] //beta of pre-exist. end. cells (Gy^-2)
#define PAR_BETA_TUM_VES   m_parD[22] //beta of neo-crea. end. cells (Gy^-2)
#define PAR_DOSE_THRES     m_parD[23] //d threshold for apoptosis (Gy)
#define PAR_ARREST_TIME    m_parD[24] //radiation-induced arrest time (h)
#define PAR_M              m_parD[25] //OER maximum value
#define PAR_K              m_parD[26] //pO2 for OER = (m + 1) / 2 (mmHg)
#define PAR_HYP_NEC_THRES  m_parD[27] //pO2 hypoxic necrosis threshold (mmHg)


class Cell: public Model{
public :
    Cell(const int i, const int j, const int l, const bool tumGrowth,
         const double doubTime, std::vector<double> cycDur, const bool res,
         const double fibDoubTime, const bool ang, const double angTime,
         const double vegfThres, std::vector<double> alpha,
         std::vector<double> beta, const double doseThres,
         const double arrestTime, const int oxy, const double hypNecThres,
         Model *const parent = 0);
    virtual ~Cell();
    virtual int calcModelOut();
    virtual int initModel();
    virtual int terminateModel();
    virtual int updateModel(const double currentTime, const double DT);
    void addToEdge(Cell *const cell);
    void calcFibProlif(double DT);
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
    std::vector<Cell *> *getEdge() const;
    bool getG0() const;
    bool getG1() const;
    bool getG2() const;
    bool getHypNec() const;
    bool getM() const;
    bool getMitCat() const;
    bool getNormVes() const;
    int getOutState() const;
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
    void setInFib(const bool input);
    void setInApop(const bool input);
    void setInHypNec(const bool input);
    void setInMitCat(const bool input);
    void setInNormVes(const bool input);
    void setInPO2(const double input);
    void setInTimer(const double input);
    void setInTum(const bool input);
    void setInTumVes(const bool input);
    void setInVegf(const double input);

protected:
    int m_i, m_j, m_l, m_r;
    Treatment *m_treatment;
    std::vector<Cell *> *m_edge;
};

bool compRedge(Cell *a, Cell *b);

#endif
