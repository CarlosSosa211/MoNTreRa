/**
 * @file constOxyTissue.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.07.17
 */

#ifndef DEF_CONSTOXYTISSUE
#define DEF_CONSTOXYTISSUE

#include <string>
#include <vector>

#include "absOxyTissue.hpp"
#include "constOxyCell.hpp"
#include "model.hpp"

//Inputs
#define CONSTOXYTISSUE_NUM_IN_B 0
#define CONSTOXYTISSUE_NUM_IN_I 0
#define CONSTOXYTISSUE_NUM_IN_D 0

//State variables
#define CONSTOXYTISSUE_NUM_ST_B 0
#define CONSTOXYTISSUE_NUM_ST_I 0
#define CONSTOXYTISSUE_NUM_ST_D 4

//Outputs
#define CONSTOXYTISSUE_NUM_OUT_B 0
#define CONSTOXYTISSUE_NUM_OUT_I 0
#define CONSTOXYTISSUE_NUM_OUT_D 5

//Internal parameters
#define CONSTOXYTISSUE_NUM_PAR_B 1
#define CONSTOXYTISSUE_NUM_PAR_I 1
#define CONSTOXYTISSUE_NUM_PAR_D 1

#define PAR_PO2_NOT_VES m_parD[0]

class ConstOxyTissue : public AbsOxyTissue{
public:
    ConstOxyTissue(const int nrow, const int ncol, const int nlayer,
                   const std::vector<bool> &inVes,
                   const std::vector<double> &inPO2, const int oxy,
                   const double hypThres, const double pO2NotVes = 0.0);
    virtual ~ConstOxyTissue();
    virtual int initModel();
    virtual int calcModelOut();
    virtual int updateModel(double currentTime, const double DT);

protected:
    int m_ncol, m_nlayer, m_nrow;
};

#endif
