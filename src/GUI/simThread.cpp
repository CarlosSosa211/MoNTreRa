#include <fstream>
#include <iostream>
#include <vector>

#include "simThread.hpp"

SimThread::SimThread(QObject *parent) : QThread(parent){
}


void SimThread::run(){
    std::ifstream fTissueDim("../OutputFilesGUI/tissueDim.dat");
    int nrow, ncol, nlayer;

    fTissueDim >> nrow >> ncol >> nlayer;

    std::cout << "Tissue dimensions: "  << std::endl;
    std::cout << nrow << " row";
    if(nrow > 1){
        std::cout << "s";
    }
    std::cout << std::endl;
    std::cout << ncol << " column";
    if(ncol > 1){
        std::cout << "s";
    }
    std::cout << std::endl;
    std::cout << nlayer << " layer";
    if(nlayer > 1){
        std::cout << "s";
    }
    std::cout << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    std::vector<bool> inTum;
    std::ifstream fInTum("inTum.dat");
    bool temp;

    fInTum >> temp;
    while(!fInTum.eof()){
        inTum.push_back(temp);
        fInTum >> temp;
    }

    fInTum.close();

    std::vector<bool> inVes;
    std::ifstream fInVes("inVes.dat");

    fInVes >> temp;
    while(!fInVes.eof()){
        inVes.push_back(temp);
        fInVes >> temp;
    }

    fInVes.close();

    std::ifstream fParam("../OutputFilesGUI/param.dat");

    std::vector<double> cycDistrib(4);

    for(int i(0); i < 4; i++){
        fParam >> cycDistrib.at(i);
    }

    bool tumGrowth;
    int edgeSize;
    double doubTime;
    std::vector<double> cycDur(4);

    fParam >> tumGrowth;
    if(tumGrowth){
        fParam >> doubTime >> edgeSize;
        for(int i(0); i < 4; i++){
            fParam >> cycDur.at(i);
        }
    }

    bool res;
    double fibDoubTime;

    fParam >> res;
    if(res){
        fParam >> fibDoubTime;
    }

    bool ang;
    double angTime, Dvegf, VmaxVegf, KmVegf;
    double hypVegf, vegfThres;

    fParam >> ang;
    if(ang){
        fParam >> angTime >> Dvegf >> VmaxVegf >> KmVegf;
        fParam >> hypVegf >> vegfThres;
    }

    bool RT;
    int schedule;
    double fraction, interval, totalDose;
    double doseThres(4.0), arrestTime;
    std::vector<double> alpha(8);
    std::vector<double> beta(8);
    Treatment *treatment;

    fParam >> RT;
    if(RT){
        for(int i(0); i < 8; i++){
            fParam >> alpha.at(i);
        }
        for(int i(0); i < 8; i++){
            fParam >> beta.at(i);
        }

        fParam >> arrestTime;

        fParam >> fraction >> totalDose >> interval >> schedule;

        treatment = new Treatment(fraction, totalDose,
                                  interval, schedule);
    }
    else{
        treatment = 0;
    }

    std::cout << treatment << std::endl;

    bool oxy;
    double constpO2, constpO2NormVes, constpO2TumVes;
    double DO2, hypThres, hypNecThres;
    double Km, pO2NormVes, pO2TumVes, Vmax;

    fParam >> oxy;
    if(oxy){
        fParam >> DO2 >> Vmax >> Km >> pO2NormVes;
        fParam >> pO2TumVes >> hypThres >> hypNecThres;
    }
    else{
        fParam >> constpO2 >> constpO2NormVes >> constpO2TumVes;
    }

    fParam.close();

    std::ifstream fSimParam("../OutputFilesGUI/simParam.dat");

    int oxySimTimeStep, simTime, simTimeStep, simType;

    fSimParam >> simType >> simTime >> simTimeStep >> oxySimTimeStep;

    fSimParam.close();

    if(simType / 2){
        double sclFac;
        Coupler *coupler;
        RootSimulator *sim;
        Tissue *model1;

        model1 = new Tissue(nrow, ncol, nlayer, inTum,
                            inVes, tumGrowth, doubTime, edgeSize,
                            cycDur, cycDistrib, res, fibDoubTime,
                            ang, angTime, vegfThres, alpha, beta,
                            doseThres, arrestTime, treatment,
                            hypNecThres);

        if(oxy){
            OxyTissue *model2;

            model2 = new OxyTissue(nrow, ncol, nlayer, inVes,
                                   Dvegf, DO2, Vmax, Km,
                                   pO2NormVes, pO2TumVes,
                                   hypThres, VmaxVegf, KmVegf,
                                   hypVegf);
            coupler = new Coupler(model1, model2);
            sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
        }

        else{
            ConstOxyTissue *model2;

            model2 = new ConstOxyTissue(nrow, ncol, nlayer,
                                        inVes, constpO2, constpO2NormVes,
                                        constpO2TumVes, hypThres);
            coupler = new Coupler(model1, model2);
            sclFac = 1.0;
        }

        sim = new RootSimulator(coupler, simTimeStep,
                                oxySimTimeStep, sclFac);

        int numIter(simTime / simTimeStep);
        double currentTime(0.0);

        sim->initSim();

        std::ofstream fTumDens("../OutputFilesGUI/tumDens.res");
        std::ofstream fKilledCells("../OutputFilesGUI/killedCells.res");
        std::ofstream fVascDens("../OutputFilesGUI/vascDens.res");
        std::ofstream fHypDens("../OutputFilesGUI/hypDens.res");
        std::ofstream fPO2Stat("../OutputFilesGUI/pO2Stat.res");
        std::ofstream fVEGFStat("../OutputFilesGUI/vegfStat.res");
        std::ofstream fState("../OutputFilesGUI/state.res");
        std::ofstream fCycle("../OutputFilesGUI/cycle.res");
        std::ofstream fTimer("../OutputFilesGUI/timer.res");
        std::ofstream fPO2("../OutputFilesGUI/po2.res");
        std::ofstream fVEGF("../OutputFilesGUI/vegf.res");

        fTumDens << currentTime << " " <<
                    model1->getOut()->at(0) << std::endl;
        fVascDens << currentTime << " " <<
                     model1->getOut()->at(6) << " " <<
                     model1->getOut()->at(7) << " " <<
                     model1->getOut()->at(8) << std::endl;
        fKilledCells << currentTime << " " <<
                        model1->getOut()->at(21) << std::endl;
        fHypDens << currentTime << " " <<
                    coupler->getModel2()->getOut()->at(0) << std::endl;
        fPO2Stat << currentTime << " " <<
                    coupler->getModel2()->getOut()->at(1) << " " <<
                    coupler->getModel2()->getOut()->at(2) << std::endl;
        fVEGFStat << currentTime << " " <<
                     coupler->getModel2()->getOut()->at(3) << " " <<
                     coupler->getModel2()->getOut()->at(4) << std::endl;
        fCycle << currentTime << " " <<
                  model1->getOut()->at(1) << " " <<
                  model1->getOut()->at(2) << " " <<
                  model1->getOut()->at(3) << " " <<
                  model1->getOut()->at(4) << " " <<
                  model1->getOut()->at(5) << std::endl;

        for(int i(0); i < model1->getNumComp(); i++){
            fState << model1->getComp()->at(i)->getOut()->at(0) << "\t";
            fTimer << model1->getComp()->at(i)->getOut()->at(1) << "\t";
            fPO2 << coupler->getModel2()->getComp()->at(i)->getOut()->at(0) << "\t";
            fVEGF << coupler->getModel2()->getComp()->at(i)->getOut()->at(1) << "\t";

        }
        fState << std::endl;
        fTimer << std::endl;
        fPO2 << std::endl;
        fVEGF << std::endl;

        emit progressMax(numIter);
        for(int j(0); j < numIter; j++){
            currentTime += simTimeStep;
            sim->simulate(currentTime, simTimeStep);

            fTumDens << currentTime << " " <<
                        model1->getOut()->at(0) << std::endl;
            fVascDens << currentTime << " " <<
                         model1->getOut()->at(6) << " " <<
                         model1->getOut()->at(7) << " " <<
                         model1->getOut()->at(8) << std::endl;
            fKilledCells << currentTime << " " <<
                            model1->getOut()->at(21) << std::endl;
            fHypDens << currentTime << " " <<
                        coupler->getModel2()->getOut()->at(0) << std::endl;
            fPO2Stat << currentTime << " " <<
                        coupler->getModel2()->getOut()->at(1) << " " <<
                        coupler->getModel2()->getOut()->at(2) << std::endl;
            fVEGFStat << currentTime << " " <<
                         coupler->getModel2()->getOut()->at(3) << " " <<
                         coupler->getModel2()->getOut()->at(4) << std::endl;
            fCycle << currentTime << " " <<
                      model1->getOut()->at(1) << " " <<
                      model1->getOut()->at(2) << " " <<
                      model1->getOut()->at(3) << " " <<
                      model1->getOut()->at(4) << " " <<
                      model1->getOut()->at(5) << std::endl;

            for(int i(0); i < model1->getNumComp(); i++){
                fState << model1->getComp()->at(i)->getOut()->at(0) << "\t";
                fTimer << model1->getComp()->at(i)->getOut()->at(1) << "\t";
                fPO2 << coupler->getModel2()->getComp()->at(i)->getOut()->at(0) << "\t";
                fVEGF << coupler->getModel2()->getComp()->at(i)->getOut()->at(1) << "\t";
            }
            fState << std::endl;
            fTimer << std::endl;
            fPO2 << std::endl;
            fVEGF << std::endl;

            emit progress(j);
        }

        fTumDens.close();
        fVascDens.close();
        fKilledCells.close();
        fHypDens.close();
        fPO2Stat.close();
        fVEGFStat.close();
        fState.close();
        fTimer.close();
        fCycle.close();
        fPO2.close();
        fVEGF.close();

        emit progress(numIter);

        sim->stop();

        std::ofstream fIntTumDens("../OutputFilesGUI/intTumDens.res");

        fIntTumDens << currentTime << " " <<
                       model1->getOut()->at(22) << std::endl;
        fIntTumDens.close();

        std::ofstream fDoseAndTime ("../OutputFilesGUI/doseAndTime.res");

        double perc[6] = {50.0, 80.0, 90.0, 95.0, 99.0, 99.9};

        for(int i(0); i < 6; i++){
            fDoseAndTime << perc[i] << " " << model1->getOut()->at(9 + i)
                         << " " << model1->getOut()->at(15 + i)
                         << std::endl;
        }

        fDoseAndTime.close();

        delete treatment;
        delete model1;
        delete coupler->getModel2();
        delete coupler;
        delete sim;
    }

    else{
        OxyTissue *model1;
        Simulator *sim;

        model1 = new OxyTissue(nrow, ncol, nlayer,
                               inVes, Dvegf, DO2, Vmax,
                               Km, pO2NormVes, pO2TumVes,
                               hypThres, VmaxVegf, KmVegf,
                               hypVegf);

        sim = new Simulator(model1, oxySimTimeStep);

        int numIter(simTime / oxySimTimeStep);
        double currentTime(0.0);

        sim->initSim();

        std::ofstream fHypDens("../OutputFilesGUI/hypDens.res");
        std::ofstream fPO2Stat("../OutputFilesGUI/pO2Stat.res");
        std::ofstream fVEGFStat("../OutputFilesGUI/vegfStat.res");
        std::ofstream fPO2("../OutputFilesGUI/po2.res");
        std::ofstream fVEGF("../OutputFilesGUI/vegf.res");

        fHypDens << currentTime << " " <<
                    model1->getOut()->at(0) << std::endl;
        fPO2Stat << currentTime << " " <<
                    model1->getOut()->at(1) << " " <<
                    model1->getOut()->at(2) << std::endl;
        fVEGFStat << currentTime << " " <<
                     model1->getOut()->at(3) << " " <<
                     model1->getOut()->at(4) << std::endl;

        for(int i(0); i < model1->getNumComp(); i++){
            fPO2 << model1->getComp()->at(i)->getOut()->at(0) << "\t";
            fVEGF << model1->getComp()->at(i)->getOut()->at(1) << "\t";
        }
        fPO2 << std::endl;
        fVEGF << std::endl;

        emit progressMax(numIter);
        for(int j(0); j < numIter; j++){
            currentTime += oxySimTimeStep;
            sim->simulate(currentTime, oxySimTimeStep);

            fHypDens << currentTime << " " <<
                        model1->getOut()->at(0) << std::endl;
            fPO2Stat << currentTime << " " <<
                        model1->getOut()->at(1) << " " <<
                        model1->getOut()->at(2) << std::endl;
            fVEGFStat << currentTime << " " <<
                         model1->getOut()->at(3) << " " <<
                         model1->getOut()->at(4) << std::endl;

            for(int i(0); i < model1->getNumComp(); i++){
                fPO2 << model1->getComp()->at(i)->getOut()->at(0) << "\t";
                fVEGF << model1->getComp()->at(i)->getOut()->at(1) << "\t";
            }
            fPO2 << std::endl;
            fVEGF << std::endl;

            emit progress(j);
        }

        emit progress(numIter);

        fHypDens.close();
        fPO2Stat.close();
        fVEGFStat.close();
        fPO2.close();
        fVEGF.close();

        sim->stop();

        delete model1;
        delete sim;
    }
    emit resultReady(simType);
}
