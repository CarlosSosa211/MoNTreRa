#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "constOxyTissue.hpp"
#include "oxyTissue.hpp"
#include "rootSimulator.hpp"
#include "sensAn.hpp"
#include "simpcell.hpp"
#include "tcp.hpp"
#include "tissue.hpp"
#include "treatment.hpp"

using namespace std;


int main(){
    //const int N(100);
    const int kp(5), L(5), p(20), P(1), N(10);
    //const int nMethod(0), nModel(1);
    //string nFRefParInt("../InputFiles/refParIntOneAlphaBeta.dat");
    //string nFRefParInt("../InputFiles/refParIntRT.dat");
    //string nFRefParInt("../InputFiles/refParIntToy.dat");
    //string nRefParMean("../InputFiles/refParMeanRT.dat");
    string nFMostRelPar("../InputFiles/mostRelParRes.dat");
    string nFLeastRelPar("../InputFiles/leastRelParRes.dat");
    string nFVarPar("../InputFiles/varParRes.dat");
    string nFInTissueDim("../InputFiles/tissueDim.dat");
    string nFInTum("../InputFiles/inTum.dat");
    string nFInVes("../InputFiles/inVes.dat");
    //vector<string> nFPar;
    //nFPar.push_back("../InputFiles/par37_2.dat");
    //nFPar.push_back("../InputFiles/par20_3.dat");
    /*string nFInTissueTCP("../InputFiles/inTissueTCP0.dat");
    string nFParTCP("../InputFiles/parTCP.dat");
    vector<string> nFTreatmentTCP;
    nFTreatmentTCP.push_back("../InputFiles/1MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/2MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/3MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/4MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/5MonFri.dat");*/

    srand(time(NULL));
    //evalR(nMethod, nModel, nFInTissueDim, nFInTum, nFInVes);
    //morrisRT(N, p, nFRefParInt, nFInTissueDim, nFInTum, nFInVes);
    //morrisToy(N, p, nFRefParInt);
    //morrisVarRangeRT(kp, L, N, p, nFRefParInt, nFInTissueDim, nFInTum, nFInVes);
    //morrisVarRangeToy(kp, L, N, p, nFRefParInt);
    //sobolRT(N, nFRefParInt, nFInTissueDim, nFInTum, nFInVes));
    //sobolToy(N, nFRefParInt);
    //sobolFromFiles(2);
    //tcp(N, nFInTissueTCP, nFParTCP, nFTreatmentTCP);
    /*tcp(N, nFInTissueTCP, nFParTCP, nFTreatmentTCP, nFInTissueDim, nFInTum,
    nFInVes);*/
    //var1ParRange(kp, L, nFRefParInt, nFInTissueDim, nFInTum, nFInVes);
    varErr(nFVarPar, nFMostRelPar, nFLeastRelPar, nFInTissueDim, nFInTum,
           nFInVes, L, P);
    //varParFromFiles(nFPar, nFInTissueDim, nFInTum, nFInVes);
    //varStoch(N, P, nFRefParInt, nFInTissueDim, nFInTum, nFInVes);
}


void model(const double *x, double *y, const int nrow,
           const int ncol, const int nlayer, const double cellSize,
           const vector<bool> &inTum, const vector<bool> & inVes){
    //int k(0);
    /*int nrow(90), ncol(90), nlayer(1);
    cout << "Creating tissue with: "  << endl;
    cout << nrow << " row";
    if(nrow > 1){
    cout << "s";
    }
    cout << endl;
    cout << ncol << " column";
    if(ncol > 1){
    cout << "s";
    }
    cout << endl;
    cout << nlayer << " layer";
    if(nlayer > 1){
    cout << "s";
    }
    cout << endl;

    cout << "---------------------------------------------" << endl;*/

    /*double tumDensRef(0.4), sigmaTumRef(0.1);
    double vascDensRef(0.03), sigmaVascRef(0.1);

    double tumDens   = tumDensRef;
    double sigmaTum  = sigmaTumRef;
    double vascDens  = vascDensRef;
    double sigmaVasc = sigmaVascRef;

    cout << "tumDens: " << tumDens << endl;
    cout << "sigmaTum: " << sigmaTum << endl;
    cout << "vascDens: " << vascDens << endl;
    cout << "sigmaVasc: " << sigmaVasc << endl;*/

    /*createInFiles(nrow, ncol, nlayer, tumDens, sigmaTum,
    vascDens, sigmaVasc, inTum, inVes);*/

    /*vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    const double tumGrowth(1.0);
    const int edgeOrder(1);
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};
    const double res(1.0);
    const double ang(1.0);

    const double tumTime(x[k]);
    k++;
    const double fibTime(x[k]);
    k++;
    const double vascTumTime(x[k]);
    k++;
    const double vegfThres(x[k]);
    k++;
    vector<double> alpha(8), beta(8);
    for(int i(0); i < 8; i++){
        alpha[i] = x[k];
        k++;
        beta[i]  = alpha[i] / x[k];
        k++;
    }*/
    /*for(int i(0); i < 8; i++){
        alpha[i] = x[k];
        beta[i]  = alpha[i] / x[k + 1];
    }
    k = k + 2;*/
    /*const double doseThres(x[k]);
    k++;
    const double arrestTime(x[k]);
    k++;
    const double hypNecThres(x[k]);
    k++;
    const double dose(x[k]);
    k++;
    double DO2(x[k]);
    k++;
    double Vmax(x[k]);
    k++;
    const double Km(x[k]);
    k++;
    double Dvegf(x[k]);
    k++;
    double VmaxVegf(x[k]);
    k++;
    const double KmVegf(x[k]);
    k++;
    const double pO2NormVes(x[k]);
    k++;
    const double pO2TumVes(x[k]);
    k++;
    const double hypThres(x[k]);
    k++;
    const double hypVegf(x[k]);
    k++;

    cout << "tumTime: "     << tumTime     << " h" <<endl;
    cout << "fibTime: "     << fibTime     << " h" << endl;
    cout << "vascTumTime: " << vascTumTime << " h" << endl;
    cout << "vegfThres: "   << vegfThres   << " mol/um^3" << endl;
    cout << "alpha: ";
    for(int i(0); i < 8; i++){
        cout << alpha[i] << " ";
    }
    cout << "Gy^-1" << endl;
    cout << "beta: ";
    for(int i(0); i < 8; i++){
        cout << beta[i] << " ";
    }
    cout << "Gy^-2" << endl;
    cout << "doseThres: "   << doseThres   << " Gy" << endl;
    cout << "arrestTime: "  << arrestTime  << " h" << endl;
    cout << "hypNecThres: " << hypNecThres << " mmHg" << endl;
    cout << "dose: "        << dose        << " Gy" << endl;
    cout << "DO2: "         << DO2         << " um^2/ms" << endl;
    cout << "Vmax: "        << Vmax        << " mmHg/ms" << endl;
    cout << "Km: "          << Km          << " mmHg" << endl;
    cout << "Dvegf: "       << Dvegf       << " um^2/ms" << endl;
    cout << "VmaxVegf: "    << VmaxVegf    << " mol/um^3ms" << endl;
    cout << "KmVegf: "      << KmVegf      << " mol/um^3" << endl;
    cout << "pO2NormVes: "  << pO2NormVes  << " mmHg" << endl;
    cout << "pO2TumVes: "   << pO2TumVes   << " mmHg" << endl;
    cout << "hypThres: "    << hypThres    << " mmHg" << endl;
    cout << "hypVegf: "     << hypVegf     << " mol/um^3" << endl;*/

    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);
    const double tumGrowth(x[k] > 0.5);
    k++;
    const double tumTime(x[k]);
    k++;
    const int edgeOrder(x[k]);
    k++;
    const double res(x[k] > 0.5);
    k++;
    const double fibTime(x[k]);
    k++;
    const double ang(x[k] > 0.5);
    k++;
    const double vascTumTime(x[k]);
    k++;
    double Dvegf(x[k]);
    k++;
    double VmaxVegf(x[k]);
    k++;
    const double KmVegf(x[k]);
    k++;
    const double vegfThres(x[k]);
    k++;
    const double hypVegf(x[k]);
    k++;
    vector<double> alpha(8), beta(8);
    for(int i(0); i < 8; i++){
        alpha[i] = x[k];
        k++;
        beta[i]  = alpha[i] / x[k];
        k++;
    }
    const double dose(x[k]);
    k++;
    const double doseThres(x[k]);
    k++;
    const double arrestTime(x[k]);
    k++;
    double oxy(x[k]);
    k++;
    const double hypNecThres(x[k]);
    k++;
    double DO2(x[k]);
    k++;
    double VmaxO2(x[k]);
    k++;
    const double KmO2(x[k]);
    k++;
    const double pO2NormVes(x[k]);
    k++;
    const double pO2TumVes(x[k]);
    k++;
    const double hypThres(x[k]);
    k++;

    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "res: "         << res         << endl;
    cout << "fibTime: "     << fibTime     << " h" << endl;
    cout << "ang: "         << ang         << endl;
    cout << "vascTumTime: " << vascTumTime << " h" << endl;
    cout << "Dvegf: "       << Dvegf       << " um^2/ms" << endl;
    cout << "VmaxVegf: "    << VmaxVegf    << " mol/um^3ms" << endl;
    cout << "KmVegf: "      << KmVegf      << " mol/um^3" << endl;
    cout << "vegfThres: "   << vegfThres   << " mol/um^3" << endl;
    cout << "hypVegf: "     << hypVegf     << " mol/um^3" << endl;
    cout << "alpha: ";
    for(int i(0); i < 8; i++){
        cout << alpha[i] << " ";
    }
    cout << "Gy^-1" << endl;
    cout << "beta: ";
    for(int i(0); i < 8; i++){
        cout << beta[i] << " ";
    }
    cout << "Gy^-2" << endl;
    cout << "dose: "        << dose        << " Gy" << endl;
    cout << "doseThres: "   << doseThres   << " Gy" << endl;
    cout << "arrestTime: "  << arrestTime  << " h" << endl;
    cout << "oxy: "         << oxy         << endl;
    cout << "hypNecThres: " << hypNecThres << " mmHg" << endl;
    cout << "DO2: "         << DO2         << " um^2/ms" << endl;
    cout << "VmaxO2: "      << VmaxO2      << " mmHg/ms" << endl;
    cout << "KmO2: "        << KmO2        << " mmHg" << endl;
    cout << "pO2NormVes: "  << pO2NormVes  << " mmHg" << endl;
    cout << "pO2TumVes: "   << pO2TumVes   << " mmHg" << endl;
    cout << "hypThres: "    << hypThres    << " mmHg" << endl;

    Treatment *treatment;
    treatment = new Treatment(dose, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, tumGrowth,
                        tumTime, edgeOrder, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, DO2, VmaxO2, KmO2,
                           pO2NormVes, pO2TumVes, hypThres);

    coupler = new Coupler(model1, model2);

    const double simTimeStep(6.0);
    const double sclFac(3.6e6 * simTimeStep / oxySimTimeStep);
    const double simTime(2160.0);
    RootSimulator *sim;

    sim = new RootSimulator(coupler, simTimeStep, oxySimTimeStep, sclFac);
    sim->initSim();
    sim->simulate(simTimeStep, simTime);
    sim->stop();

    double endTreatTumDens(model1->getOut()->at(24));
    double threeMonTumDens(model1->getOut()->at(25));
    double tumVol(model1->getOut()->at(23));
    double intTumDens(model1->getOut()->at(22));
    double killed50(model1->getOut()->at(28));
    double killed80(model1->getOut()->at(29));
    double killed90(model1->getOut()->at(30));
    double killed95(model1->getOut()->at(31));
    double timeTo95(model1->getOut()->at(12));
    double killed99(model1->getOut()->at(32));
    double timeTo99(model1->getOut()->at(13));
    double killed999(model1->getOut()->at(33));
    double rec(model1->getOut()->at(34));
    double recTumDens(model1->getOut()->at(26));
    double recTime(model1->getOut()->at(27));

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "endTreatTumDens: " << endTreatTumDens << endl;
    cout << "3MonTumDens: "     << threeMonTumDens << endl;
    cout << "tumVol: "          << tumVol << " mm3" << endl;
    cout << "intTumDens: "      << intTumDens << endl;
    cout << "killed50: "        << killed50 << endl;
    cout << "killed80: "        << killed80 << endl;
    cout << "killed90: "        << killed90 << endl;
    cout << "killed95: "        << killed95 << endl;
    cout << "timeTo95: "        << timeTo95 << " h" << endl;
    cout << "killed99: "        << killed99 << endl;
    cout << "timeTo99: "        << timeTo99 << " h" << endl;
    cout << "killed999: "       << killed999 << endl;
    cout << "rec: "             << rec << endl;
    cout << "recTumDens: "      << recTumDens << endl;
    cout << "recTime: "         << recTime << " h" << endl;

    y[0]  = endTreatTumDens;
    y[1]  = threeMonTumDens;
    y[2]  = tumVol;
    y[3]  = intTumDens;
    y[4]  = killed50;
    y[5]  = killed80;
    y[6]  = killed90;
    y[7]  = killed95;
    y[8]  = timeTo95;
    y[9]  = killed99;
    y[10] = timeTo99;
    y[11] = killed999;
    y[12] = rec;
    y[13] = recTumDens;
    y[14] = recTime;
}


void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize, const vector<bool> & inVes,
           const string nFPO2){
    int k(0);
    const double ang(x[k] > 0.5);
    k++;
    double Dvegf(x[k]);
    k++;
    double VmaxVegf(x[k]);
    k++;
    const double KmVegf(x[k]);
    k++;
    const double hypVegf(x[k]);
    k++;
    double oxy(x[k]);
    k++;
    double DO2(x[k]);
    k++;
    double VmaxO2(x[k]);
    k++;
    const double KmO2(x[k]);
    k++;
    const double pO2NormVes(x[k]);
    k++;
    const double pO2TumVes(x[k]);
    k++;
    const double hypThres(x[k]);
    k++;

    cout << "ang: "         << ang         << endl;
    cout << "Dvegf: "       << Dvegf       << " um^2/ms" << endl;
    cout << "VmaxVegf: "    << VmaxVegf    << " mol/um^3ms" << endl;
    cout << "KmVegf: "      << KmVegf      << " mol/um^3" << endl;
    cout << "hypVegf: "     << hypVegf     << " mol/um^3" << endl;
    cout << "oxy: "         << oxy         << endl;
    cout << "DO2: "         << DO2         << " um^2/ms" << endl;
    cout << "VmaxO2: "      << VmaxO2      << " mmHg/ms" << endl;
    cout << "KmO2: "        << KmO2        << " mmHg" << endl;
    cout << "pO2NormVes: "  << pO2NormVes  << " mmHg" << endl;
    cout << "pO2TumVes: "   << pO2TumVes   << " mmHg" << endl;
    cout << "hypThres: "    << hypThres    << " mmHg" << endl;

    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

    OxyTissue *model1;
    Simulator *sim;

    model1 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, DO2, VmaxO2, KmO2,
                           pO2NormVes, pO2TumVes, hypThres);
    sim = new Simulator(model1, oxySimTimeStep);

    const double simTime(1000.0);

    sim->initSim();
    sim->simulate(oxySimTimeStep, simTime);
    sim->stop();

    std::ofstream fPO2(nFPO2.c_str());

    for(int i(0); i < model1->getNumComp(); i++){
        fPO2 << model1->getComp()->at(i)->getOut()->at(0) << "\t";
    }

    delete model1;
    delete sim;
}


void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize, const vector<bool> &inTum,
           const vector<bool> & inVes, const string nFTumDens,
           const string nFTumVol, const string nFVascDens,
           const string nFKilledCells, const string nFDeadDens,
           const string nFCycle, const string nFHypDens, const string nFPO2Stat,
           const string nFVegfStat){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);
    const double tumGrowth(x[k]);
    k++;
    const double tumTime(x[k]);
    k++;
    const int edgeOrder(x[k]);
    k++;
    const double res(x[k]);
    k++;
    const double fibTime(x[k]);
    k++;
    const double ang(x[k]);
    k++;
    const double vascTumTime(x[k]);
    k++;
    double Dvegf(x[k]);
    k++;
    double VmaxVegf(x[k]);
    k++;
    const double KmVegf(x[k]);
    k++;
    const double vegfThres(x[k]);
    k++;
    const double hypVegf(x[k]);
    k++;
    vector<double> alpha(8), beta(8);
    for(int i(0); i < 8; i++){
        alpha[i] = x[k];
        k++;
        beta[i]  = alpha[i] / x[k];
        k++;
    }
    const double dose(x[k]);
    k++;
    const double doseThres(x[k]);
    k++;
    const double arrestTime(x[k]);
    k++;
    double oxy(x[k]);
    k++;
    const double hypNecThres(x[k]);
    k++;
    double DO2(x[k]);
    k++;
    double VmaxO2(x[k]);
    k++;
    const double KmO2(x[k]);
    k++;
    const double pO2NormVes(x[k]);
    k++;
    const double pO2TumVes(x[k]);
    k++;
    const double hypThres(x[k]);
    k++;

    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "res: "         << res         << endl;
    cout << "fibTime: "     << fibTime     << " h" << endl;
    cout << "ang: "         << ang         << endl;
    cout << "vascTumTime: " << vascTumTime << " h" << endl;
    cout << "Dvegf: "       << Dvegf       << " um^2/ms" << endl;
    cout << "VmaxVegf: "    << VmaxVegf    << " mol/um^3ms" << endl;
    cout << "KmVegf: "      << KmVegf      << " mol/um^3" << endl;
    cout << "vegfThres: "   << vegfThres   << " mol/um^3" << endl;
    cout << "hypVegf: "     << hypVegf     << " mol/um^3" << endl;
    cout << "alpha: ";
    for(int i(0); i < 8; i++){
        cout << alpha[i] << " ";
    }
    cout << "Gy^-1" << endl;
    cout << "beta: ";
    for(int i(0); i < 8; i++){
        cout << beta[i] << " ";
    }
    cout << "Gy^-2" << endl;
    cout << "dose: "        << dose        << " Gy" << endl;
    cout << "doseThres: "   << doseThres   << " Gy" << endl;
    cout << "arrestTime: "  << arrestTime  << " h" << endl;
    cout << "oxy: "         << oxy         << endl;
    cout << "hypNecThres: " << hypNecThres << " mmHg" << endl;
    cout << "DO2: "         << DO2         << " um^2/ms" << endl;
    cout << "VmaxO2: "      << VmaxO2      << " mmHg/ms" << endl;
    cout << "KmO2: "        << KmO2        << " mmHg" << endl;
    cout << "pO2NormVes: "  << pO2NormVes  << " mmHg" << endl;
    cout << "pO2TumVes: "   << pO2TumVes   << " mmHg" << endl;
    cout << "hypThres: "    << hypThres    << " mmHg" << endl;

    Treatment *treatment;
    treatment = new Treatment(dose, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, tumGrowth,
                        tumTime, edgeOrder, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double simTimeStep(6.0);
    const double oxySimTimeStep(10.0);
    double sclFac;

    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;
    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;

    if(oxy == 0){
        ConstOxyTissue *model2;

        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, 0.0, 0.0, 0.0,
                                    0.0);
        coupler = new Coupler(model1, model2);
        sclFac = 1.0;
    }

    else if(oxy == 1){
        OxyTissue *model2;

        model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                               VmaxVegf, KmVegf, hypVegf, DO2, VmaxO2, KmO2,
                               pO2NormVes, pO2TumVes, hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
    }

    const double simTime(2160.0);
    RootSimulator *sim;

    sim = new RootSimulator(coupler, simTimeStep,
                            oxySimTimeStep, sclFac);

    double currentTime(0.0);
    ofstream fTumDens(nFTumDens), fTumVol(nFTumVol);
    ofstream fKilledCells(nFKilledCells), fVascDens(nFVascDens);
    ofstream fDeadDens(nFDeadDens), fCycle(nFCycle), fHypDens(nFHypDens);
    ofstream  fPO2Stat(nFPO2Stat), fVEGFStat(nFVegfStat);

    sim->initSim();

    fTumDens     << currentTime << " " << model1->getOut()->at(0) << endl;
    fTumVol      << currentTime << " " << model1->getOut()->at(23) << endl;
    fVascDens    << currentTime << " " << model1->getOut()->at(6) << " " <<
                    model1->getOut()->at(7) << " " <<
                    model1->getOut()->at(8) << endl;
    fKilledCells << currentTime << " " << model1->getOut()->at(21) << endl;
    fDeadDens    << currentTime << " " << model1->getOut()->at(37) << endl;
    fHypDens     << currentTime << " " <<
                    coupler->getModel2()->getOut()->at(0) << endl;
    fPO2Stat     << currentTime << " " <<
                    coupler->getModel2()->getOut()->at(1) << " " <<
                    coupler->getModel2()->getOut()->at(2) << endl;
    fVEGFStat    << currentTime << " " <<
                    coupler->getModel2()->getOut()->at(3) << " " <<
                    coupler->getModel2()->getOut()->at(4) << endl;
    fCycle       << currentTime << " " << model1->getOut()->at(1) << " " <<
                    model1->getOut()->at(2) << " " <<
                    model1->getOut()->at(3) << " " <<
                    model1->getOut()->at(4) << " " <<
                    model1->getOut()->at(5) << endl;

    int numIter(simTime / simTimeStep);

    for(int j(0); j < numIter; j++){
        currentTime += simTimeStep;
        sim->simulate(currentTime, simTimeStep);

        fTumDens     << currentTime << " " << model1->getOut()->at(0) << endl;
        fTumVol      << currentTime << " " << model1->getOut()->at(23) << endl;
        fVascDens    << currentTime << " " << model1->getOut()->at(6) << " " <<
                        model1->getOut()->at(7) << " " <<
                        model1->getOut()->at(8) << endl;
        fKilledCells << currentTime << " " << model1->getOut()->at(21) << endl;
        fDeadDens    << currentTime << " " << model1->getOut()->at(37) << endl;
        fHypDens     << currentTime << " " <<
                        coupler->getModel2()->getOut()->at(0) << endl;
        fPO2Stat     << currentTime << " " <<
                        coupler->getModel2()->getOut()->at(1) << " " <<
                        coupler->getModel2()->getOut()->at(2) << endl;
        fVEGFStat    << currentTime << " " <<
                        coupler->getModel2()->getOut()->at(3) << " " <<
                        coupler->getModel2()->getOut()->at(4) << endl;
        fCycle       << currentTime << " " << model1->getOut()->at(1) << " " <<
                        model1->getOut()->at(2) << " " <<
                        model1->getOut()->at(3) << " " <<
                        model1->getOut()->at(4) << " " <<
                        model1->getOut()->at(5) << endl;
    }

    fTumDens.close();
    fTumVol.close();
    fVascDens.close();
    fKilledCells.close();
    fDeadDens.close();
    fHypDens.close();
    fPO2Stat.close();
    fVEGFStat.close();

    sim->stop();

    double endTreatTumDens(model1->getOut()->at(24));
    double threeMonTumDens(model1->getOut()->at(25));
    double tumVol(model1->getOut()->at(23));
    double intTumDens(model1->getOut()->at(22));
    double killed50(model1->getOut()->at(28));
    double killed80(model1->getOut()->at(29));
    double killed90(model1->getOut()->at(30));
    double killed95(model1->getOut()->at(31));
    double timeTo95(model1->getOut()->at(12));
    double killed99(model1->getOut()->at(32));
    double timeTo99(model1->getOut()->at(13));
    double killed999(model1->getOut()->at(33));
    double rec(model1->getOut()->at(34));
    double recTumDens(model1->getOut()->at(26));
    double recTime(model1->getOut()->at(27));

    delete treatment;
    delete model1;
    delete coupler->getModel2();
    delete coupler;
    delete sim;

    cout << "endTreatTumDens: " << endTreatTumDens << endl;
    cout << "3MonTumDens: "     << threeMonTumDens << endl;
    cout << "tumVol: "          << tumVol << " mm3" << endl;
    cout << "intTumDens: "      << intTumDens << endl;
    cout << "killed50: "        << killed50 << endl;
    cout << "killed80: "        << killed80 << endl;
    cout << "killed90: "        << killed90 << endl;
    cout << "killed95: "        << killed95 << endl;
    cout << "timeTo95: "        << timeTo95 << " h" << endl;
    cout << "killed99: "        << killed99 << endl;
    cout << "timeTo99: "        << timeTo99 << " h" << endl;
    cout << "killed999: "       << killed999 << endl;
    cout << "rec: "             << rec << endl;
    cout << "recTumDens: "      << recTumDens << endl;
    cout << "recTime: "         << recTime << " h" << endl;

    y[0]  = endTreatTumDens;
    y[1]  = threeMonTumDens;
    y[2]  = tumVol;
    y[3]  = intTumDens;
    y[4]  = killed50;
    y[5]  = killed80;
    y[6]  = killed90;
    y[7]  = killed95;
    y[8]  = timeTo95;
    y[9]  = killed99;
    y[10] = timeTo99;
    y[11] = killed999;
    y[12] = rec;
    y[13] = recTumDens;
    y[14] = recTime;
}


void oxy(const int N, const string nFInTissueOxy, const string nFParOxy,
         const string nFInTissueDim, const string nFInVes){
    const int K(12), nOut(1);
    bool art(0);
    int nrow, ncol, nlayer;
    double cellSize, vascDens, sigmaVasc;
    vector<bool> inVes;

    readInFilesOxy(nFInTissueOxy, art, nrow, ncol, nlayer, cellSize, vascDens,
                   sigmaVasc);

    if(!art){
        readInFiles(nFInTissueDim, nFInVes, nrow, ncol, nlayer, cellSize,
                    inVes);
    }

    double x[K], y[nOut];
    ifstream fParOxy(nFParOxy.c_str());

    for(int k(0); k < K; k++){
        fParOxy >> x[k];
    }
    fParOxy.close();

    for(int j(0); j < N; j++){
        string nFPO2("../OutputFiles/pO2_" + to_string(j) + ".res");
        if(art){
            createInFiles(nrow, ncol, nlayer, vascDens, sigmaVasc, inVes);
        }
        model(x, y, nrow, ncol, nlayer, cellSize, inVes, nFPO2);
        cout << j + 1 << " out of " << N << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
    }
}


void toyModel(double *x, double *y){
    y[0] = x[0] + 2.0 * x[1] + x[2] * x[2] + x[3] * x[4];
}

