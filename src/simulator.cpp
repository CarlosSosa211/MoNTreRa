/**
 * @file simulator.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include "simulator.hpp"

using namespace std;

Simulator::Simulator(){
    m_model = 0;
    m_currentTime = 0.0;
    m_DT = 1; //h
}


Simulator::Simulator(Model *model, const double DT){
    m_model = model;
    m_currentTime = 0.0;
    m_DT = DT;
}


Simulator::~Simulator(){
}


void Simulator::setModel(Model *model){
    m_model = model;
}


void Simulator::initSim(){ 
    m_model->initModel();
    m_model->calcModelOut();
}


void Simulator::simulate(const double currentTime,
                         const double simTime){
    int numIter(simTime / m_DT);

    m_currentTime = currentTime;
    for(int j(0); j < numIter; j++){
        if(m_model->updateModel(m_currentTime, m_DT)){
            break;
        }
        m_model->calcModelOut();
        m_currentTime += m_DT;
    }
}


void Simulator::stop(){
    m_model->terminateModel();
    m_model->calcModelOut();
}
