/**
 * @file rootSimulator.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef DEF_ROOTSIMULATOR
#define DEF_ROOTSIMULATOR

#include "coupler.hpp"
#include "model.hpp"
#include "simulator.hpp"

class RootSimulator{
public:
    RootSimulator();
    RootSimulator(Coupler *coupler, const double DT1, const double DT2,
                  const double sclFac);
    ~RootSimulator();
    void initSim();
    void simulate(const double currentTime, const double simTime);
    void stop();
private:
    double m_DT1, m_DT2, m_sclFac;
    double m_currentTime;
    Coupler *m_coupler;
    Model *m_model1, *m_model2;
    Simulator *m_sim1, *m_sim2;
};

#endif
