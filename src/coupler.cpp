/**
 * @file coupler.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include "coupler.hpp"

using namespace std;

/*------------------------------------------------------------------------------
 * Constructor of the class Coupler.
 *
 * Inputs:
 *  - model1: pointer to model 1,
 *  - model2: pointer to model 2.
------------------------------------------------------------------------------*/

Coupler::Coupler(Model *model1, Model *model2) : Model(0, 0, 0, 0, 0, 0, 0, 0,
                                                       0, 0, 0, 0, 2){
    m_comp->at(0) = model1;
    m_comp->at(1) = model2;
}


/*------------------------------------------------------------------------------
 * Destructor of the class Coupler.
------------------------------------------------------------------------------*/

Coupler::~Coupler(){
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model calcModelOut method.
------------------------------------------------------------------------------*/

int Coupler::calcModelOut(){
    return 0;
}

/*------------------------------------------------------------------------------
 * Redefinition of the Model initModel method.
------------------------------------------------------------------------------*/

int Coupler::initModel(){
    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model terminalModel method.
------------------------------------------------------------------------------*/

int Coupler::terminateModel(){
    return 0;
}


/*------------------------------------------------------------------------------
 * Redefinition of the Model updateModel method.
 *
 * Inputs:
 *  - currentTime: simulation current time,
 *  - DT: simulation timestep.
------------------------------------------------------------------------------*/

int Coupler::updateModel(const double currentTime, const double DT){
    if(getModel1()->getNumComp() ==  getModel2()->getNumComp()){
        for(int k(0); k < getModel1()->getNumComp(); k++){
            ((Cell *)getModel1()->getComp()->at(k))->
                    setInPO2(getModel2()->getComp()->at(k)->getOutD()[0]);
            ((Cell *)getModel1()->getComp()->at(k))->
                    setInVegf(getModel2()->getComp()->at(k)->getOutD()[1]);

            /*((Cell *)getModel1()->getComp()->at(k))->
                    setInPO2(getModel2()->getOutD()[2]);
            ((Cell *)getModel1()->getComp()->at(k))->
                    setInVegf(getModel2()->getOutD()[4]);*/

            ((AbsOxyCell *)getModel2()->getComp()->at(k))->
                    setInNormVes(((Cell *)getModel1()->getComp()->at(k))->
                                 getNormVes());

            ((AbsOxyCell *)getModel2()->getComp()->at(k))->
                    setInTumVes(((Cell *)getModel1()->getComp()->at(k))->
                                getTumVes());

            ((AbsOxyCell *)getModel2()->getComp()->at(k))->
                    setInDead(((Cell *)getModel1()->getComp()->at(k))->
                              getDead());
        }
    }
    return 0;
}


/*------------------------------------------------------------------------------
 * This function gets the model 1.
------------------------------------------------------------------------------*/

Model *Coupler::getModel1() const{
    return  m_comp->at(0);
}


/*------------------------------------------------------------------------------
 * This function gets the model 2.
------------------------------------------------------------------------------*/

Model *Coupler::getModel2() const{
    return  m_comp->at(1);
}
