/**
 * @file simulator.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#ifndef DEF_SIMULATOR
#define DEF_SIMULATOR

#include "model.hpp"

class Simulator{	
public:
  Simulator();
  Simulator(Model *model, const double DT);
  ~Simulator();
  void initSim();
  void setModel(Model *model);
  void simulate(const double currentTime, const double simTime);
  void stop();
  
private:
  double m_currentTime, m_DT;
  Model *m_model;
};

#endif
