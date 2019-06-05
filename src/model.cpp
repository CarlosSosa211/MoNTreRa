/**
 * @file model.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include "model.hpp"
#include <iostream>

using namespace std;

Model::Model(const int numInB, const int numInI, const int numInD,
             const int numStB, const int numStI, const int numStD,
             const int numOutB, const int numOutI, const int numOutD,
             const int numParB, const int numParI, const int numParD,
             const int numComp){
    m_numIn[0] = numInB + numInI + numInD;
    m_numIn[1] = numInB;
    m_numIn[2] = numInI;
    m_numIn[3] = numInD;
    m_numSt[0] = numStB + numStI + numStD;
    m_numSt[1] = numStB;
    m_numSt[2] = numStI;
    m_numSt[3] = numStD;
    m_numOut[0] = numOutB + numOutI + numOutD;
    m_numOut[1] = numOutB;
    m_numOut[2] = numOutI;
    m_numOut[3] = numOutD;
    m_numPar[0] = numParB + numParI + numParD;
    m_numPar[1] = numParB;
    m_numPar[2] = numParI;
    m_numPar[3] = numParD;
    m_numComp = numComp;

    m_inB  = new bool[numInB];
    m_inI  = new int[numInI];
    m_inD  = new double[numInD];
    m_stB  = new bool[numStB];
    m_stI  = new int[numStI];
    m_stD  = new double[numStD];
    m_outB = new bool[numOutB];
    m_outI = new int[numOutI];
    m_outD = new double[numOutD];
    m_parB = new bool[numParB];
    m_parI = new int[numParI];
    m_parD = new double[numParD];
    m_comp = new vector<Model *>((unsigned int)numComp, 0);
}


Model::~Model(){
    delete [] m_inB;
    delete [] m_inI;
    delete [] m_inD;
    delete [] m_stB;
    delete [] m_stI;
    delete [] m_stD;
    delete [] m_outB;
    delete [] m_outI;
    delete [] m_outD;
    delete [] m_parB;
    delete [] m_parI;
    delete [] m_parD;
    delete m_comp;
}



vector<Model*> *Model::getComp() const{
    return m_comp;
}


bool *Model::getInB() const{
    return m_inB;
}


double *Model::getInD() const{
    return m_inD;
}


int *Model::getInI() const{
    return m_inI;
}


bool *Model::getOutB() const{
    return m_outB;
}


double *Model::getOutD() const{
    return m_outD;
}


int *Model::getOutI() const{
    return m_outI;
}


bool *Model::getParB() const{
    return m_parB;
}


int *Model::getParI() const{
    return m_parI;
}


double *Model::getParD() const{
    return m_parD;
}


bool *Model::getStB() const{
    return m_stB;
}


int *Model::getStI() const{
    return m_stI;
}


double *Model::getStD() const{
    return m_stD;
}


int Model::getNumIn() const{
    return m_numIn[0];
}


int Model::getNumSt() const{
    return m_numSt[0];
}


int Model::getNumOut() const{
    return m_numOut[0];

}


int Model::getNumPar() const{
    return m_numPar[0];
}


int Model::getNumComp() const{
    return m_numComp;
}
