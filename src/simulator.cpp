/**
 * @file simulator.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.19.17
 */

#include "simulator.hpp"

using namespace std;

/*------------------------------------------------------------------------------
 * Constructor of the class Simulator.
 *
 * Inputs:
 *  - model: pointer to the model to be simulated,
 *  - DT: simulation timestep.
------------------------------------------------------------------------------*/

Simulator::Simulator(Model *model, const double DT){
    m_model = model;
    m_currentTime = 0.0;
    m_DT = DT;
}


/*------------------------------------------------------------------------------
 * Destructor of the class Simulator.
------------------------------------------------------------------------------*/

Simulator::~Simulator(){
}


/*------------------------------------------------------------------------------
 * This function sets the model.
 *
 * Inputs:
 *  - model: pointer to the model to be simulated.
------------------------------------------------------------------------------*/

void Simulator::setModel(Model *model){
    m_model = model;
}


/*------------------------------------------------------------------------------
 * This function initialises the simulation.
------------------------------------------------------------------------------*/
void Simulator::initSim(){ 
    m_model->initModel();
    m_model->calcModelOut();
}


/*------------------------------------------------------------------------------
 * This function runs a simulation of the model.
 *
 * Inputs:
 *  - currentTime: initial simulation time,
 *  - simTime: duration of the simulation.
------------------------------------------------------------------------------*/

void Simulator::simulate(const double currentTime, const double simTime){
    int numIter(simTime / m_DT);

    m_currentTime = currentTime;
    for(int j(0); j < numIter; j++){
        if(m_model->updateModel(m_currentTime, m_DT)){
            j = numIter;
        }
        m_model->calcModelOut();
        m_currentTime += m_DT;
    }
}


/*------------------------------------------------------------------------------
 * This function stops the simulation.
------------------------------------------------------------------------------*/

void Simulator::stop(){
    m_model->terminateModel();
    m_model->calcModelOut();
}
