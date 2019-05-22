
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

/*------------------------------------------------------------------------------
 * Constructor of the class RootSimulator.
 *
 * Inputs:
 *  - coupler: pointer to the coupler containing the 2 models to be simulated,
 *  - DT1: simulation 1 timestep,
 *  - DT2: simulation 2 timestep,
 *  - sclFac: scale factor between the 2 timesteps.
------------------------------------------------------------------------------*/

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


/*------------------------------------------------------------------------------
 * Destructor of the class RootSimulator.
------------------------------------------------------------------------------*/

RootSimulator::~RootSimulator(){
    delete m_sim1;
    delete m_sim2;
}


/*------------------------------------------------------------------------------
 * This function initialises the simulation.
------------------------------------------------------------------------------*/

void RootSimulator::initSim(){
    m_sim1->initSim();
    m_sim2->initSim();
}


/*------------------------------------------------------------------------------
 * This function runs a simulation of the models.
 *
 * Inputs:
 *  - currentTime: initial simulation time,
 *  - simTime: duration of the simulation.
------------------------------------------------------------------------------*/

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


/*------------------------------------------------------------------------------
 * This function stops the simulation.
------------------------------------------------------------------------------*/

void RootSimulator::stop(){
    m_sim1->stop();
    m_sim2->stop();
}
