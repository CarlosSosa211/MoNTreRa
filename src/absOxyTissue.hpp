/**
 * @file absOxyTissue.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.22.19
 */

#ifndef DEF_ABSOXYTISSUE
#define DEF_ABSOXYTISSUE

#include <string>
#include <vector>

#include "absOxyCell.hpp"
#include "model.hpp"

//State variables
#define ST_OXYTISSUE_VES_DENS      m_stD[0] //vascular density
#define ST_OXYTISSUE_NORM_VES_DENS m_stD[1] //pre-existing vascular density
#define ST_OXYTISSUE_TUM_VES_DENS  m_stD[2] //neo-created vascular density
#define ST_OXYTISSUE_DEAD_DENS     m_stD[3] //dead cell density

//Outputs
#define OUT_HYP_DENS  m_outD[0] //hypoxic density
#define OUT_PO2_MED   m_outD[1] //median pO2 (mmHg)
#define OUT_PO2_MEAN  m_outD[2] //mean pO2 (mmHg)
#define OUT_VEGF_MED  m_outD[3] //median VEGF concentration (mol/um^3)
#define OUT_VEGF_MEAN m_outD[4] //mean VEGF concentration (mol/um^3)

//Internal parameters
#define PAR_OXYTISSUE_ANG m_parB[0] //angiogenesis

#define PAR_OXYTISSUE_OXY m_parI[0] //oxygenation

class AbsOxyTissue : public Model{
public:
    AbsOxyTissue(const int numInB, const int numInI, const int numInD,
                 const int numStB, const int numStI, const int numStD,
                 const int numOutB, const int numOutI, const int numOutD,
                 const int numParB, const int numParI, const int numParD,
                 const int numComp);
    virtual ~AbsOxyTissue();
    virtual int initModel() = 0;
    virtual int calcModelOut() = 0;
    virtual int terminateModel();
    virtual int updateModel(double currentTime, const double DT) = 0;
    int getNumDead() const;
    int getNumHyp() const;
    int getNumNormVes() const;
    int getNumTumVes() const;
    int getNumVes() const;

protected:
    int m_ncol, m_nlayer, m_nrow;
};

#endif
