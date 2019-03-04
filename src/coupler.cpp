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

Coupler::Coupler(Model *model1, Model *model2) : Model(0, 0, 1, 0, 2){
    m_comp->at(0) = model1;
    m_comp->at(1) = model2;
}


Coupler::~Coupler(){
}


int Coupler::calcModelOut(){
    return 0;
}


int Coupler::initModel(){
    return 0;
}


int Coupler::terminateModel(){
    return 0;
}


int Coupler::updateModel(const double currentTime,
                         const double DT){
    if(getModel1()->getNumComp() ==  getModel2()->getNumComp()){
        for(int k(0); k < getModel1()->getNumComp(); k++){
            ((Cell *)getModel1()->getComp()->at(k))->
                    setInPO2(getModel2()->getComp()->at(k)->getOut()->at(0));

            ((Cell *)getModel1()->getComp()->at(k))->
                    setInVegf(getModel2()->getComp()->at(k)->getOut()->at(1));

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


Model *Coupler::getModel1() const{
    return  m_comp->at(0);
}


Model *Coupler::getModel2() const{
    return  m_comp->at(1);
}
