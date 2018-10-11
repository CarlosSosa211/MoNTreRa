/**
 * @file coupler.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef DEF_COUPLER
#define DEF_COUPLER

#include "model.hpp"
#include "oxyTissue.hpp"
#include "tissue.hpp"

class Coupler : public Model{
public :
    Coupler(Model *model1, Model *model2);
    virtual ~Coupler();
    virtual int calcModelOut();
    virtual int initModel();
    virtual int terminateModel();
    virtual int updateModel(const double currentTime = 0,
                            const double DT = 0);
    Model *getModel1() const;
    Model *getModel2() const;
};

#endif
