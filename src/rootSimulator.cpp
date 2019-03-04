
/**
 * @file rootSimulator.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include <fstream>

#include "rootSimulator.hpp"

using namespace std;

RootSimulator::RootSimulator(){
}


RootSimulator::RootSimulator(Coupler *coupler, const double DT1,
                             const double DT2, const double sclFac){
    m_coupler = coupler;
    m_model1 = m_coupler->getModel1();
    m_model2 = m_coupler->getModel2();
    m_DT1 = DT1;
    m_DT2 = DT2;
    m_sclFac = sclFac;
    m_currentTime = 0.0;
    m_sim1 = new Simulator(m_model1, DT1);
    m_sim2 = new Simulator(m_model2, DT2);
}


RootSimulator::~RootSimulator(){
    delete m_sim1;
    delete m_sim2;
}


void RootSimulator::initSim(){
    m_sim1->initSim();
    m_sim2->initSim();
}


void RootSimulator::simulate(const double currentTime, const double simTime){
    int numIter(simTime / m_DT1);
    m_currentTime = currentTime;

    for(int j(0); j < numIter; j++){
        m_sim2->simulate(m_currentTime, m_sclFac);
        m_coupler->updateModel();
        m_sim1->simulate(m_currentTime, m_DT1);
        m_coupler->updateModel();
        m_currentTime += m_DT1;
    }
}


void RootSimulator::stop(){
    m_sim1->stop();
    m_sim2->stop();
}
