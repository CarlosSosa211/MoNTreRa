/**
 * @file absOxyCell.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.07.17
 */

#ifndef DEF_ABSTOXYCELL
#define DEF_ABSTOXYCELL

#include "model.hpp"

//Inputs
#define IN_OXYCELL_DEAD     m_inB[0] //dead cell
#define IN_OXYCELL_NORM_VES m_inB[1] //pre-existing endothelial cell
#define IN_OXYCELL_TUM_VES  m_inB[2] //neo-created endothelial cell

//State variables
#define ST_OXYCELL_DEAD        m_stB[0] //dead cell
#define ST_OXYCELL_NORM_VES    m_stB[1] //pre-existing endothelial cell
#define ST_OXYCELL_TUM_VES     m_stB[2] //neo-created endothelial cell
#define ST_OXYCELL_VES         m_stB[3] //endothelial
#define ST_HYP                 m_stB[4] //hypoxic cell
#define ST_OXYCELL_OXY_STABLE  m_stB[5] //pO2 stability
#define ST_OXYCELL_VEGF_STABLE m_stB[6] //VEGF concentration stability

#define ST_OXYCELL_PO2  m_stD[0] //pO2 (mmHg)
#define ST_OXYCELL_VEGF m_stD[1] //VEGF concentration (mol/um^3)

//Outputs
#define OUT_PO2  m_outD[0] //pO2 (mmHg)
#define OUT_VEGF m_outD[1] //VEGF concentration (mol/um^3)

//Internal parameters
#define PAR_HYP_VEGF     m_parD[0] //fixed VEGF conc. for hyp. cells (mol/um^3)
#define PAR_PO2_NORM_VES m_parD[1] //fixed pO2 val. of pre-ex. end. cells (mmHg)
#define PAR_PO2_TUM_VES  m_parD[2] //fixed pO2 val. of neo-cr. end. cells (mmHg)
#define PAR_HYP_THRES    m_parD[3] //pO2 hypoxia threshold (mmHg)


class AbsOxyCell: public Model{
public :
    AbsOxyCell(const int numInB, const int numInI, const int numInD,
               const int numStB, const int numStI, const int numStD,
               const int numOutB, const int numOutI, const int numOutD,
               const int numParB, const int numParI, const int numParD);
    virtual ~AbsOxyCell();
    virtual int calcModelOut();
    virtual int initModel();
    virtual int terminateModel();
    virtual int updateModel(const double currentTime, const double DT) = 0;
    bool getDead() const;
    bool getHyp() const;
    double getOutPO2() const;
    double getOutVEGF() const;
    bool getOxyStable() const;
    double getPO2() const;
    double getVEGF() const;
    bool getVegfStable () const;
    bool getNormVes() const;
    bool getTumVes() const;
    bool getVes() const;
    void setInDead(const bool input);
    void setInNormVes(const bool input);
    void setInTumVes(const bool input);
};

#endif
