/**
 * @file constOxyTissue.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.07.17
 */

#ifndef DEF_CONSTOXYTISSUE
#define DEF_CONSTOXYTISSUE

#include <string>
#include <vector>

#include "constOxyCell.hpp"
#include "model.hpp"

//Outputs
#define OUT_HYP_DENS  m_out->at(0)
#define OUT_PO2_MED   m_out->at(1)
#define OUT_PO2_MEAN  m_out->at(2)
#define OUT_VEGF_MED  m_out->at(3)
#define OUT_VEGF_MEAN m_out->at(4)

class ConstOxyTissue : public Model{
public:
    ConstOxyTissue(const int nrow, const int ncol, const int nlayer,
                   const std::string nFInVes, const double pO2,
                   const double pO2NormVes, const double pO2TumVes,
                   const double hypThres);
    ConstOxyTissue(const int nrow, const int ncol, const int nlayer,
                   const std::vector<bool> &inVes, const double pO2,
                   const double pO2NormVes, const double pO2TumVes,
                   const double hypThres);
    virtual ~ConstOxyTissue();
    virtual int initModel();
    virtual int calcModelOut();
    virtual int updateModel(double currentTime, const double DT);
    virtual int terminateModel();
    int getNumHyp() const;

protected:
    int m_ncol, m_nlayer, m_nrow;
};

#endif
