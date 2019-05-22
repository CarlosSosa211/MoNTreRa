/**
 * @file model.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef DEF_MODEL
#define DEF_MODEL

#include <vector>

class Model {
public :
    Model(const int numInB, const int numInI, const int numInD,
          const int numStB, const int numStI, const int numStD,
          const int numOutB, const int numOutI, const int numOutD,
          const int numParB, const int numParI, const int numParD,
          const int numComp);
    virtual ~Model();
    virtual int calcModelOut() = 0;
    virtual int initModel() = 0;
    virtual int terminateModel() = 0;
    virtual int updateModel(const double currentTime = 0, const double DT = 0) =
            0;
    bool *getInB() const;
    bool *getStB() const;
    bool *getOutB() const;
    bool *getParB() const;
    int *getInI() const;
    int *getStI() const;
    int *getOutI() const;
    int *getParI() const;
    double *getInD() const;
    double *getStD() const;
    double *getOutD() const;
    double *getParD() const;
    std::vector<Model*> *getComp() const;
    int getNumIn() const;
    int getNumSt() const;
    int getNumOut() const;
    int getNumPar() const;
    int getNumComp() const;

protected:
    int m_numIn[4], m_numSt[4], m_numOut[4], m_numPar[4], m_numComp;
    bool *m_inB, *m_stB, *m_outB, *m_parB;
    int *m_inI, *m_stI, *m_outI, *m_parI;
    double *m_inD, *m_stD, *m_outD, *m_parD;
    std::vector<Model *> *m_comp;
    Model *m_parent;
};

#endif
