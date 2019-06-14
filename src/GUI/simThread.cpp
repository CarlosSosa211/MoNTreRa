#include <fstream>
#include <iostream>
#include <vector>

#include "simThread.hpp"

/*------------------------------------------------------------------------------
 * Constructor of the class SimThread.
------------------------------------------------------------------------------*/

SimThread::SimThread(QObject *parent) : QThread(parent){
}


/*------------------------------------------------------------------------------
 * Redefinition of the QThread run method.
------------------------------------------------------------------------------*/

void SimThread::run(){
    std::ifstream fTissueDim("../OutputFilesGUI/tissueDim.dat");
    int nrow, ncol, nlayer;
    double cellSize;

    fTissueDim >> nrow >> ncol >> nlayer >> cellSize;

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
    int edgeOrder;
    double doubTime;
    std::vector<double> cycDur(4);

    fParam >> tumGrowth;
    if(tumGrowth){
        fParam >> doubTime >> edgeOrder;
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
        fParam >> vegfThres >> hypVegf;
    }

    bool RT;
    int schedule;
    double fraction, interval, totalDose;
    double doseThres, arrestTime;
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

        fParam >> fraction >> totalDose >> interval >> schedule;
        fParam >> doseThres >> arrestTime;

        treatment = new Treatment(fraction, totalDose, interval, schedule);
    }
    else{
        treatment = 0;
    }

    std::cout << treatment << std::endl;

    int oxy;
    double DO2, VmaxO2, KmO2;
    double hypThres, hypNecThres;
    double pO2NormVes, pO2TumVes;
    double constpO2NotVes, constpO2NormVes, constpO2TumVes;

    fParam >> oxy;
    switch(oxy){
    case 1:{
        fParam >> hypNecThres >> DO2 >> VmaxO2 >> KmO2;
        fParam >> pO2NormVes >> pO2TumVes >> hypThres;
        break;
    }

    case 2:{
        fParam >> hypNecThres >> DO2 >> VmaxO2 >> KmO2;
        fParam >> pO2NormVes >> pO2TumVes >> hypThres;
        break;
    }

    case 3:{
        fParam >> constpO2NormVes >> constpO2TumVes >> hypThres;
        break;
    }

    case 4:{
        fParam >> constpO2NotVes >> constpO2NormVes >> constpO2TumVes >>
                hypThres;
        break;
    }
    }

    fParam.close();

    int oxySimTimeStep, simTime, simTimeStep, simType;
    std::ifstream fSimParam("../OutputFilesGUI/simParam.dat");

    fSimParam >> simType >> simTime >> simTimeStep >> oxySimTimeStep;

    fSimParam.close();

    if(simType / 2){
        double sclFac;
        Coupler *coupler;
        RootSimulator *sim;
        Tissue *model1;
        AbsOxyTissue *model2;

        model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes,
                            tumGrowth, doubTime, edgeOrder, cycDur, cycDistrib,
                            res, fibDoubTime, ang, angTime, vegfThres, alpha,
                            beta, treatment, doseThres, arrestTime, oxy,
                            hypNecThres);
        switch(oxy){
        case 0:{
            model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy);
            coupler = new Coupler(model1, model2);
            sclFac = oxySimTimeStep;
            break;
        }

        case 1:{
            model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang,
                                   Dvegf, VmaxVegf, KmVegf, hypVegf, oxy, DO2,
                                   VmaxO2, KmO2, pO2NormVes, pO2TumVes,
                                   hypThres);
            coupler = new Coupler(model1, model2);
            sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
            break;
        }

        case 2:{
            model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang,
                                   Dvegf, VmaxVegf, KmVegf, hypVegf, oxy, DO2,
                                   VmaxO2, KmO2, pO2NormVes, pO2TumVes,
                                   hypThres);
            coupler = new Coupler(model1, model2);
            sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
            break;
        }

        case 3:{
            model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy,
                                        constpO2NormVes, constpO2TumVes,
                                        hypThres);
            coupler = new Coupler(model1, model2);
            sclFac = oxySimTimeStep;
            break;
        }

        case 4:{
            model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy,
                                        constpO2NormVes, constpO2TumVes,
                                        hypThres, constpO2NotVes);
            coupler = new Coupler(model1, model2);
            sclFac = oxySimTimeStep;
            break;
        }
        }

        sim = new RootSimulator(coupler, simTimeStep, oxySimTimeStep, sclFac);

        int numIter(simTime / simTimeStep);
        double currentTime(0.0);

        sim->initSim();

        std::ofstream fTumDens("../OutputFilesGUI/tumDens.res");
        std::ofstream fTumVol("../OutputFilesGUI/tumVol.res");
        std::ofstream fKilledCells("../OutputFilesGUI/killedCells.res");
        std::ofstream fVascDens("../OutputFilesGUI/vascDens.res");
        std::ofstream fHypDens("../OutputFilesGUI/hypDens.res");
        std::ofstream fDeadDens("../OutputFilesGUI/deadCellsDens.res");
        std::ofstream fPO2Stat("../OutputFilesGUI/pO2Stat.res");
        std::ofstream fVEGFStat("../OutputFilesGUI/vegfStat.res");
        std::ofstream fState("../OutputFilesGUI/state.res");
        std::ofstream fCycle("../OutputFilesGUI/cycle.res");
        std::ofstream fTimer("../OutputFilesGUI/timer.res");
        std::ofstream fPO2("../OutputFilesGUI/po2.res");
        std::ofstream fVEGF("../OutputFilesGUI/vegf.res");

        fTumDens     << currentTime << " " << model1->getOutD()[18] <<
                        std::endl;
        fTumVol      << currentTime << " " << model1->getOutD()[19] <<
                        std::endl;
        fVascDens    << currentTime << " " << model1->getOutD()[20] <<
                        " " << model1->getOutD()[21] <<
                        " " << model1->getOutD()[22] << std::endl;
        fKilledCells << currentTime << " " << model1->getOutD()[23] <<
                        std::endl;
        fHypDens     << currentTime << " " << model2->getOutD()[0] << std::endl;
        fDeadDens    << currentTime << " " << model1->getOutD()[24] <<
                        std::endl;
        fPO2Stat     << currentTime << " " <<
                        coupler->getModel2()->getOutD()[1] << " " <<
                        coupler->getModel2()->getOutD()[2] << std::endl;
        fVEGFStat    << currentTime << " " <<
                        coupler->getModel2()->getOutD()[3] << " " <<
                        coupler->getModel2()->getOutD()[4] << std::endl;
        fCycle       << currentTime << " " << model1->getOutD()[25] << " " <<
                        model1->getOutD()[26] << " " <<
                        model1->getOutD()[27] << " " <<
                        model1->getOutD()[28] << " " <<
                        model1->getOutD()[29] << std::endl;

        for(int i(0); i < model1->getNumComp(); i++){
            fState << model1->getComp()->at(i)->getOutI()[0] << "\t";
            fTimer << model1->getComp()->at(i)->getOutI()[1] << "\t";
            fPO2   << model2->getComp()->at(i)->getOutD()[0] << "\t";
            fVEGF  << model2->getComp()->at(i)->getOutD()[1] << "\t";
        }
        fState << std::endl;
        fTimer << std::endl;
        fPO2   << std::endl;
        fVEGF  << std::endl;

        emit progressMax(numIter);
        for(int j(0); j < numIter; j++){
            currentTime += simTimeStep;
            sim->simulate(currentTime, simTimeStep);

            fTumDens     << currentTime << " " << model1->getOutD()[18] <<
                            std::endl;
            fTumVol      << currentTime << " " << model1->getOutD()[19] <<
                            std::endl;
            fVascDens    << currentTime << " " << model1->getOutD()[20] <<
                            " " << model1->getOutD()[21] <<
                            " " << model1->getOutD()[22] << std::endl;
            fKilledCells << currentTime << " " << model1->getOutD()[23] <<
                            std::endl;
            fHypDens     << currentTime << " " << model2->getOutD()[0] <<
                            std::endl;
            fDeadDens    << currentTime << " " << model1->getOutD()[24] <<
                            std::endl;
            fPO2Stat     << currentTime << " " << model2->getOutD()[1] << " " <<
                            model2->getOutD()[2] << std::endl;
            fVEGFStat    << currentTime << " " << model2->getOutD()[3] << " " <<
                            model2->getOutD()[4] << std::endl;
            fCycle       << currentTime << " " << model1->getOutD()[25] <<
                            " " << model1->getOutD()[26] <<
                            " " << model1->getOutD()[27] <<
                            " " << model1->getOutD()[28] <<
                            " " << model1->getOutD()[29] << std::endl;

            for(int i(0); i < model1->getNumComp(); i++){
                fState << model1->getComp()->at(i)->getOutI()[0] << "\t";
                fTimer << model1->getComp()->at(i)->getOutI()[1] << "\t";
                fPO2   << model2->getComp()->at(i)->getOutD()[0] << "\t";
                fVEGF  << model2->getComp()->at(i)->getOutD()[1] << "\t";

            }
            fState << std::endl;
            fTimer << std::endl;
            fPO2   << std::endl;
            fVEGF  << std::endl;

            emit progress(j);
        }

        fTumDens.close();
        fTumVol.close();
        fVascDens.close();
        fKilledCells.close();
        fHypDens.close();
        fDeadDens.close();
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

        fIntTumDens << currentTime << " " << model1->getOutD()[2] << std::endl;
        fIntTumDens.close();

        std::ofstream fPercKilled("../OutputFilesGUI/percKilled.res");

        double perc[6] = {50.0, 80.0, 90.0, 95.0, 99.0, 99.9};

        for(int i(0); i < 6; i++){
            fPercKilled << perc[i] << " " << model1->getOutB()[0 + i]
                        << " " << model1->getOutD()[3 + i]
                        << " " << model1->getOutD()[9 + i] << std::endl;
        }

        fPercKilled.close();

        std::ofstream fEndTreatTumDens("../OutputFilesGUI/endTreatTumDens.res");
        fEndTreatTumDens << treatment->getDuration() << " " <<
                            model1->getOutD()[0] << std::endl;
        fEndTreatTumDens.close();

        std::ofstream f3MonTumDens("../OutputFilesGUI/3MonTumDens.res");
        f3MonTumDens << 2160.0 << " " << model1->getOutD()[1] << std::endl;
        f3MonTumDens.close();

        std::ofstream fRec("../OutputFilesGUI/rec.res");
        fRec << model1->getOutB()[6] << " " << model1->getOutD()[17] <<
                                        " " << model1->getOutD()[16] <<
                                        std::endl;
        fRec.close();

        delete treatment;
        delete model1;
        delete model2;
        delete coupler;
        delete sim;
    }

    else{
        OxyTissue *model1;
        Simulator *sim;

        model1 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                               VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2,
                               KmO2, pO2NormVes, pO2TumVes, hypThres);
        sim = new Simulator(model1, oxySimTimeStep);

        int numIter(simTime / oxySimTimeStep);
        double currentTime(0.0);

        sim->initSim();

        std::ofstream fHypDens("../OutputFilesGUI/hypDens.res");
        std::ofstream fPO2Stat("../OutputFilesGUI/pO2Stat.res");
        std::ofstream fVEGFStat("../OutputFilesGUI/vegfStat.res");
        std::ofstream fPO2("../OutputFilesGUI/po2.res");
        std::ofstream fVEGF("../OutputFilesGUI/vegf.res");

        fHypDens  << currentTime << " " << model1->getOutD()[0] << std::endl;
        fPO2Stat  << currentTime << " " << model1->getOutD()[1] << " " <<
                     model1->getOutD()[2] << std::endl;
        fVEGFStat << currentTime << " " << model1->getOutD()[3] << " " <<
                     model1->getOutD()[4] << std::endl;

        for(int i(0); i < model1->getNumComp(); i++){
            fPO2  << model1->getComp()->at(i)->getOutD()[0] << "\t";
            fVEGF << model1->getComp()->at(i)->getOutD()[1] << "\t";
        }
        fPO2  << std::endl;
        fVEGF << std::endl;

        emit progressMax(numIter);
        for(int j(0); j < numIter; j++){
            currentTime += oxySimTimeStep;
            sim->simulate(currentTime, oxySimTimeStep);

            fHypDens  << currentTime << " " << model1->getOutD()[0] <<
                         std::endl;
            fPO2Stat  << currentTime << " " << model1->getOutD()[1] << " " <<
                         model1->getOutD()[2] << std::endl;
            fVEGFStat << currentTime << " " << model1->getOutD()[3] << " " <<
                         model1->getOutD()[4] << std::endl;

            for(int i(0); i < model1->getNumComp(); i++){
                fPO2  << model1->getComp()->at(i)->getOutD()[0] << "\t";
                fVEGF << model1->getComp()->at(i)->getOutD()[1] << "\t";
            }
            fPO2  << std::endl;
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

        std::ofstream fOxyStable("../OutputFilesGUI/oxyStable.res");

        fOxyStable << model1->getOutB()[0] << " " << model1->getOutD()[5] <<
                                              std::endl;
        fOxyStable.close();

        std::ofstream fVegfStable("../OutputFilesGUI/vegfStable.res");

        fVegfStable << model1->getOutB()[1] << " " << model1->getOutD()[6] <<
                                               std::endl;
        fVegfStable.close();

        delete model1;
        delete sim;
    }
    emit resultReady(simType);
}
