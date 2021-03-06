﻿#include "evalModel.hpp"

using namespace std;

/*------------------------------------------------------------------------------
 * This functions evaluates the model of tumour growth and response to
 * radiotherapy to obtain its scalar outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize, const vector<bool> &inTum,
           const vector<bool> &inVes){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);
    const int edgeOrder(x[k]);
    k++;
    const bool tumGrowth(x[k] > 0.5);
    k++;
    const double tumTime(x[k]);
    k++;
    const bool res(x[k] > 0.5);
    k++;
    const double fibTime(x[k]);
    k++;
    const bool ang(x[k] > 0.5);
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
    //const double dose(x[k]);
    //k++;
    const double doseThres(x[k]);
    k++;
    const double arrestTime(x[k]);
    k++;
    const int oxy(x[k]);
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

    cout << "nrow: "        << nrow        << endl;
    cout << "ncol: "        << ncol        << endl;
    cout << "nlayer: "      << nlayer      << endl;
    cout << "cellSize: "    << cellSize    << endl;
    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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
    //cout << "dose: "        << dose        << " Gy" << endl;
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
    treatment = new Treatment(2.0, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2, KmO2,
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

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This functions evaluates the model of tumour growth and response to
 * radiotherapy using an artificial initial tissue to obtain its scalar outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);
    const double tumDens(x[k]);
    k++;
    const double sigmaTum(x[k]);
    k++;
    const double vascDens(x[k]);
    k++;
    const double sigmaVasc(x[k]);
    k++;
    const int edgeOrder(x[k]);
    k++;
    const bool tumGrowth(x[k] > 0.5);
    k++;
    const double tumTime(x[k]);
    k++;
    const bool res(x[k] > 0.5);
    k++;
    const double fibTime(x[k]);
    k++;
    const bool ang(x[k] > 0.5);
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
    const double doseThres(x[k]);
    k++;
    const double arrestTime(x[k]);
    k++;
    const int oxy(x[k]);
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

    cout << "tumDens: "     << tumDens     << endl;
    cout << "sigmaTum: "    << sigmaTum    << endl;
    cout << "vascDens: "    << vascDens    << endl;
    cout << "sigmaVasc: "   << sigmaVasc   << endl;
    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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
    treatment = new Treatment(2.0, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    vector<bool> inTum, inVes;
    createInFiles(nrow, ncol, nlayer, tumDens, sigmaTum, vascDens, sigmaVasc,
                  inTum, inVes);

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2, KmO2,
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

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This functions evaluates the model of oxygenation and angiogenesis to obtain
 * pO2(x, t) values.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - nFPO2: name of the file where pO2(x, t) values will be written.
------------------------------------------------------------------------------*/

void model(const double *x, const int nrow, const int ncol, const int nlayer,
           const double cellSize, const vector<bool> & inVes,
           const string nFPO2){
    int k(0);
    const bool ang(x[k] > 0.5);
    k++;
    double Dvegf(x[k]);
    k++;
    double VmaxVegf(x[k]);
    k++;
    const double KmVegf(x[k]);
    k++;
    const double hypVegf(x[k]);
    k++;
    const int oxy(x[k]);
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
                           VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2, KmO2,
                           pO2NormVes, pO2TumVes, hypThres);
    sim = new Simulator(model1, oxySimTimeStep);

    const double simTime(1000.0);

    sim->initSim();
    sim->simulate(oxySimTimeStep, simTime);
    sim->stop();

    ofstream fPO2(nFPO2.c_str());

    for(int i(0); i < model1->getNumComp(); i++){
        fPO2 << model1->getComp()->at(i)->getOutD()[0] << "\t";
    }

    delete model1;
    delete sim;
}


/*------------------------------------------------------------------------------
 * This functions evaluates the model of tumour growth and response to
 * radiotherapy to obtain both its scalar and time-dependent outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - nFTumDens: name of the file where the tumour density values will be
 *  written,
 *  - nFTumVol: name of the file where the tumour volume values will be written,
 *  - nFVascDens: name of the file where the vascular density values will be
 *  written,
 *  - nFKilledCells: name of the file where the killed tumour cell percentages
 *  values will be written,
 *  - nFDeadDens: name of the file where the dead cell density values will be
 *  written,
 *  - nFCycle: name of the file where the cycle distributon values will be
 *  written,
 *  - nFHypDens: name of the file where the hypoxic density values will be
 *  written,
 *  - nFPO2Stat: name of the file where the pO2 statistic values will be
 *  written,
 *  - nFVegfStat: name of the file where the pO2 statistic values will be
 *  written.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize, const vector<bool> &inTum,
           const vector<bool> &inVes, const string nFTumDens,
           const string nFTumVol, const string nFVascDens,
           const string nFKilledCells, const string nFDeadDens,
           const string nFCycle, const string nFHypDens, const string nFPO2Stat,
           const string nFVegfStat){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);
    const int edgeOrder(x[k]);
    k++;
    const bool tumGrowth(x[k]);
    k++;
    const double tumTime(x[k]);
    k++;
    const bool res(x[k]);
    k++;
    const double fibTime(x[k]);
    k++;
    const bool ang(x[k]);
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
    const int oxy(x[k]);
    k++;
    const double hypNecThres(x[k]);
    k++;

    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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

    Treatment *treatment;
    treatment = new Treatment(dose, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    AbsOxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double simTimeStep(6.0);
    const double oxySimTimeStep(10.0);
    double sclFac;

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;

    switch(oxy){
    case 0:{
        cout << treatment;

        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy);
        coupler = new Coupler(model1, model2);
        sclFac = 1.0;
        break;
    }

    case 1:{
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

        cout << "DO2: "        << DO2        << " um^2/ms" << endl;
        cout << "VmaxO2: "     << VmaxO2     << " mmHg/ms" << endl;
        cout << "KmO2: "       << KmO2       << " mmHg" << endl;
        cout << "pO2NormVes: " << pO2NormVes << " mmHg" << endl;
        cout << "pO2TumVes: "  << pO2TumVes  << " mmHg" << endl;
        cout << "hypThres: "   << hypThres   << " mmHg" << endl;
        cout << treatment;

        DO2      *= oxySimTimeStep;
        VmaxO2   *= oxySimTimeStep;

        model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang,
                               Dvegf, VmaxVegf, KmVegf, hypVegf, oxy, DO2,
                               VmaxO2, KmO2, pO2NormVes, pO2TumVes, hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
        break;
    }

    case 2:{
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

        cout << "DO2: "        << DO2        << " um^2/ms" << endl;
        cout << "VmaxO2: "     << VmaxO2     << " mmHg/ms" << endl;
        cout << "KmO2: "       << KmO2       << " mmHg" << endl;
        cout << "pO2NormVes: " << pO2NormVes << " mmHg" << endl;
        cout << "pO2TumVes: "  << pO2TumVes  << " mmHg" << endl;
        cout << "hypThres: "   << hypThres   << " mmHg" << endl;
        cout << treatment;

        DO2      *= oxySimTimeStep;
        VmaxO2   *= oxySimTimeStep;

        model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang,
                               Dvegf, VmaxVegf, KmVegf, hypVegf, oxy, DO2,
                               VmaxO2, KmO2, pO2NormVes, pO2TumVes, hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
        break;
    }

    case 3:{
        const double constPO2NormVes(x[k]);
        k++;
        const double constPO2TumVes(x[k]);
        k++;
        const double hypThres(x[k]);
        k++;
        cout << treatment;

        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy,
                                    constPO2NormVes, constPO2TumVes, hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = oxySimTimeStep;

        break;
    }

    case 4:{
        const double constPO2NotVes(x[k]);
        k++;
        const double constPO2NormVes(x[k]);
        k++;
        const double constPO2TumVes(x[k]);
        k++;
        const double hypThres(x[k]);
        k++;
        cout << treatment;

        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy,
                                    constPO2NormVes, constPO2TumVes, hypThres,
                                    constPO2NotVes);
        coupler = new Coupler(model1, model2);
        sclFac = oxySimTimeStep;
        break;
    }
    }

    const double simTime(3600.0);
    RootSimulator *sim;
    sim = new RootSimulator(coupler, simTimeStep, oxySimTimeStep, sclFac);
    double currentTime(0.0);
    ofstream fTumDens(nFTumDens), fTumVol(nFTumVol);
    ofstream fKilledCells(nFKilledCells), fVascDens(nFVascDens);
    ofstream fDeadDens(nFDeadDens), fCycle(nFCycle), fHypDens(nFHypDens);
    ofstream fPO2Stat(nFPO2Stat), fVEGFStat(nFVegfStat);

    sim->initSim();

    fTumDens     << currentTime << " " << model1->getOutD()[24] << endl;
    fTumVol      << currentTime << " " << model1->getOutD()[25] << endl;
    fVascDens    << currentTime << " " << model1->getOutD()[26] <<
                    " " << model1->getOutD()[27] <<
                    " " << model1->getOutD()[28] << endl;
    fKilledCells << currentTime << " " << model1->getOutD()[29] << endl;
    fHypDens     << currentTime << " " << model2->getOutD()[0] << endl;
    fDeadDens    << currentTime << " " << model1->getOutD()[30] << endl;
    fPO2Stat     << currentTime << " " << model2->getOutD()[1] << " " <<
                    model2->getOutD()[2] << endl;
    fVEGFStat    << currentTime << " " << model2->getOutD()[3] << " " <<
                    model2->getOutD()[4] << endl;
    fCycle       << currentTime << " " << model1->getOutD()[31] <<
                    " " << model1->getOutD()[32] <<
                    " " << model1->getOutD()[33] <<
                    " " << model1->getOutD()[34] <<
                    " " << model1->getOutD()[35] << endl;

    int numIter(simTime / simTimeStep);

    for(int j(0); j < numIter; j++){
        currentTime += simTimeStep;
        sim->simulate(currentTime, simTimeStep);

        fTumDens     << currentTime << " " << model1->getOutD()[24] << endl;
        fTumVol      << currentTime << " " << model1->getOutD()[25] << endl;
        fVascDens    << currentTime << " " << model1->getOutD()[26] <<
                        " " << model1->getOutD()[27] <<
                        " " << model1->getOutD()[28] << endl;
        fKilledCells << currentTime << " " << model1->getOutD()[29] << endl;
        fHypDens     << currentTime << " " << model2->getOutD()[0]  << endl;
        fDeadDens    << currentTime << " " << model1->getOutD()[30] << endl;
        fPO2Stat     << currentTime << " " << model2->getOutD()[1] << " " <<
                        model2->getOutD()[2] << endl;
        fVEGFStat    << currentTime << " " << model2->getOutD()[3] << " " <<
                        model2->getOutD()[4] << endl;
        fCycle       << currentTime << " " << model1->getOutD()[31] <<
                        " " << model1->getOutD()[32] <<
                        " " << model1->getOutD()[33] <<
                        " " << model1->getOutD()[34] <<
                        " " << model1->getOutD()[35] << endl;
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

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This functions evaluates the model of tumour growth and response to
 * radiotherapy to obtain both its scalar and time-dependent outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - treatment: pointer to the treatment to be considered,
 *  - nFTumDens: name of the file where the tumour density values will be
 *  written,
 *  - nFTumVol: name of the file where the tumour volume values will be written,
 *  - nFVascDens: name of the file where the vascular density values will be
 *  written,
 *  - nFKilledCells: name of the file where the killed tumour cell percentages
 *  values will be written,
 *  - nFDeadDens: name of the file where the dead cell density values will be
 *  written,
 *  - nFCycle: name of the file where the cycle distributon values will be
 *  written,
 *  - nFHypDens: name of the file where the hypoxic density values will be
 *  written,
 *  - nFPO2Stat: name of the file where the pO2 statistic values will be
 *  written,
 *  - nFVegfStat: name of the file where the pO2 statistic values will be
 *  written,
 *  - nFState: name of the file where the type of every cell will be written,
 *  - nFTimer: name of the file where the timer of every cell in its cycle will
 *  be written,
 *  - nFPO2: name of the file where the pO2 value of every cell will be
 *  written,
 *  - nFVegf: name of the file where the VEGF concentration value of every cell
 *  will be written.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize, const vector<bool> &inTum,
           const vector<bool> &inVes, Treatment *const treatment,
           const string nFTumDens, const string nFTumVol,
           const string nFVascDens, const string nFKilledCells,
           const string nFDeadDens, const string nFCycle,
           const string nFHypDens, const string nFPO2Stat,
           const string nFVegfStat, const string nFState, const string nFTimer,
           const string nFPO2, const string nFVegf){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);
    const int edgeOrder(x[k]);
    k++;
    const bool tumGrowth(x[k]);
    k++;
    const double tumTime(x[k]);
    k++;
    const bool res(x[k]);
    k++;
    const double fibTime(x[k]);
    k++;
    const bool ang(x[k]);
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
    const double doseThres(x[k]);
    k++;
    const double arrestTime(x[k]);
    k++;
    const int oxy(x[k]);
    k++;
    const double hypNecThres(x[k]);
    k++;

    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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
    cout << "doseThres: "   << doseThres   << " Gy" << endl;
    cout << "arrestTime: "  << arrestTime  << " h" << endl;
    cout << "oxy: "         << oxy         << endl;
    cout << "hypNecThres: " << hypNecThres << " mmHg" << endl;

    Coupler *coupler;
    Tissue *model1;
    AbsOxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double simTimeStep(6.0);
    const double oxySimTimeStep(10.0);
    double sclFac;

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;

    switch(oxy){
    case 0:{
        cout << treatment;

        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy);
        coupler = new Coupler(model1, model2);
        sclFac = 1.0;
        break;
    }

    case 1:{
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

        cout << "DO2: "        << DO2        << " um^2/ms" << endl;
        cout << "VmaxO2: "     << VmaxO2     << " mmHg/ms" << endl;
        cout << "KmO2: "       << KmO2       << " mmHg" << endl;
        cout << "pO2NormVes: " << pO2NormVes << " mmHg" << endl;
        cout << "pO2TumVes: "  << pO2TumVes  << " mmHg" << endl;
        cout << "hypThres: "   << hypThres   << " mmHg" << endl;
        cout << treatment;

        DO2      *= oxySimTimeStep;
        VmaxO2   *= oxySimTimeStep;

        model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang,
                               Dvegf, VmaxVegf, KmVegf, hypVegf, oxy, DO2,
                               VmaxO2, KmO2, pO2NormVes, pO2TumVes, hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
        break;
    }

    case 2:{
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

        cout << "DO2: "        << DO2        << " um^2/ms" << endl;
        cout << "VmaxO2: "     << VmaxO2     << " mmHg/ms" << endl;
        cout << "KmO2: "       << KmO2       << " mmHg" << endl;
        cout << "pO2NormVes: " << pO2NormVes << " mmHg" << endl;
        cout << "pO2TumVes: "  << pO2TumVes  << " mmHg" << endl;
        cout << "hypThres: "   << hypThres   << " mmHg" << endl;
        cout << treatment;

        DO2      *= oxySimTimeStep;
        VmaxO2   *= oxySimTimeStep;

        model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang,
                               Dvegf, VmaxVegf, KmVegf, hypVegf, oxy, DO2,
                               VmaxO2, KmO2, pO2NormVes, pO2TumVes, hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
        break;
    }

    case 3:{
        const double constPO2NormVes(x[k]);
        k++;
        const double constPO2TumVes(x[k]);
        k++;
        const double hypThres(x[k]);
        k++;
        cout << treatment;

        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy,
                                    constPO2NormVes, constPO2TumVes, hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = oxySimTimeStep;

        break;
    }

    case 4:{
        const double constPO2NotVes(x[k]);
        k++;
        const double constPO2NormVes(x[k]);
        k++;
        const double constPO2TumVes(x[k]);
        k++;
        const double hypThres(x[k]);
        k++;
        cout << treatment;

        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, oxy,
                                    constPO2NormVes, constPO2TumVes, hypThres,
                                    constPO2NotVes);
        coupler = new Coupler(model1, model2);
        sclFac = oxySimTimeStep;
        break;
    }
    }

    const double simTime(2160.0);
    RootSimulator *sim;
    sim = new RootSimulator(coupler, simTimeStep, oxySimTimeStep, sclFac);
    double currentTime(0.0);
    ofstream fTumDens(nFTumDens), fTumVol(nFTumVol);
    ofstream fKilledCells(nFKilledCells), fVascDens(nFVascDens);
    ofstream fDeadDens(nFDeadDens), fCycle(nFCycle), fHypDens(nFHypDens);
    ofstream fPO2Stat(nFPO2Stat), fVegfStat(nFVegfStat);
    ofstream fState(nFState), fTimer(nFTimer), fPO2(nFPO2), fVegf(nFVegf);

    sim->initSim();

    fTumDens     << currentTime << " " << model1->getOutD()[24] << endl;
    fTumVol      << currentTime << " " << model1->getOutD()[25] << endl;
    fVascDens    << currentTime << " " << model1->getOutD()[26] <<
                    " " << model1->getOutD()[27] <<
                    " " << model1->getOutD()[28] << endl;
    fKilledCells << currentTime << " " << model1->getOutD()[29] << endl;
    fHypDens     << currentTime << " " << model2->getOutD()[0] << endl;
    fDeadDens    << currentTime << " " << model1->getOutD()[30] << endl;
    fPO2Stat     << currentTime << " " << model2->getOutD()[1] << " " <<
                    model2->getOutD()[2] << endl;
    fVegfStat    << currentTime << " " << model2->getOutD()[3] << " " <<
                    model2->getOutD()[4] << endl;
    fCycle       << currentTime << " " << model1->getOutD()[31] <<
                    " " << model1->getOutD()[32] <<
                    " " << model1->getOutD()[33] <<
                    " " << model1->getOutD()[34] <<
                    " " << model1->getOutD()[35] << endl;

    for(int i(0); i < model1->getNumComp(); i++){
        fState << model1->getComp()->at(i)->getOutI()[0] << "\t";
        fTimer << model1->getComp()->at(i)->getOutI()[1] << "\t";
        fPO2   << model2->getComp()->at(i)->getOutD()[0] << "\t";
        fVegf  << model2->getComp()->at(i)->getOutD()[1] << "\t";
    }
    fState << endl;
    fTimer << endl;
    fPO2   << endl;
    fVegf  << endl;

    int numIter(simTime / simTimeStep);

    for(int j(0); j < numIter; j++){
        currentTime += simTimeStep;
        sim->simulate(currentTime, simTimeStep);

        fTumDens     << currentTime << " " << model1->getOutD()[24] << endl;
        fTumVol      << currentTime << " " << model1->getOutD()[25] << endl;
        fVascDens    << currentTime << " " << model1->getOutD()[26] <<
                        " " << model1->getOutD()[27] <<
                        " " << model1->getOutD()[28] << endl;
        fKilledCells << currentTime << " " << model1->getOutD()[29] << endl;
        fHypDens     << currentTime << " " << model2->getOutD()[0]  << endl;
        fDeadDens    << currentTime << " " << model1->getOutD()[30] << endl;
        fPO2Stat     << currentTime << " " << model2->getOutD()[1] << " " <<
                        model2->getOutD()[2] << endl;
        fVegfStat    << currentTime << " " << model2->getOutD()[3] << " " <<
                        model2->getOutD()[4] << endl;
        fCycle       << currentTime << " " << model1->getOutD()[31] <<
                        " " << model1->getOutD()[32] <<
                        " " << model1->getOutD()[33] <<
                        " " << model1->getOutD()[34] <<
                        " " << model1->getOutD()[35] << endl;

        for(int i(0); i < model1->getNumComp(); i++){
            fState << model1->getComp()->at(i)->getOutI()[0] << "\t";
            fTimer << model1->getComp()->at(i)->getOutI()[1] << "\t";
            fPO2   << model2->getComp()->at(i)->getOutD()[0] << "\t";
            fVegf  << model2->getComp()->at(i)->getOutD()[1] << "\t";

        }
        fState << endl;
        fTimer << endl;
        fPO2   << endl;
        fVegf  << endl;
    }

    fTumDens.close();
    fTumVol.close();
    fVascDens.close();
    fKilledCells.close();
    fDeadDens.close();
    fHypDens.close();
    fPO2Stat.close();
    fVegfStat.close();

    sim->stop();

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);

    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This function evaluates N times the model of oxygenation and angiogenesis to
 * obtain pO2(x, t) values, using parameters and a tissue architecture defined
 * in input files.
 *
 * Inputs:
 *  - N: number of repetitions,
 *  - nFInTissueOxy: name of the file containing the tissue architecture,
 *  - nFParOxy: name of the file containing the values of the parameters of the
 *  model,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen; it is empty if an artificial tissue is used,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration; it is empty if an artificial tissue is used.
------------------------------------------------------------------------------*/

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
        readInFiles(nFInTissueDim, nFInVes, nrow, ncol, nlayer,
                    cellSize, inVes);
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
        model(x, nrow, ncol, nlayer, cellSize, inVes, nFPO2);
        cout << j + 1 << " out of " << N << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
    }
}

/*------------------------------------------------------------------------------
 * This functions evaluates the reduced model of tumour growth and response to
 * radiotherapy to obtain its scalar outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void reducedFracModel(const double *x, double *y, const int nrow, const int ncol,
                  const int nlayer, const double cellSize,
                  const vector<bool> &inTum, const vector<bool> &inVes){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    const double tumTime(330.0);
    const int edgeOrder(1);
    const bool tumGrowth(true);
    const bool res(false);
    const double fibTime(0.0);
    const bool ang(false);
    const double vascTumTime(0.0);
    double Dvegf(0.0);
    double VmaxVegf(0.0);
    const double KmVegf(0.0);
    const double vegfThres(0.0);
    const double hypVegf(0.0);
    vector<double> alpha(8), beta(8);
    alpha[0] = 0.0;
    alpha[1] = 0.120;
    alpha[2] = 0.111;
    alpha[3] = 0.165;
    alpha[4] = 0.184;
    alpha[5] = 0.15;
    alpha[6] = 0.0;
    alpha[7] = 0.0;
    beta[0]  = 0.0;
    beta[1]  = 0.120 / 5.5;
    beta[2]  = 0.111 / 5.5;
    beta[3]  = 0.165 / 5.5;
    beta[4]  = 0.184 / 5.5;
    beta[5] = 0.15 / 5.5;
    beta[6] = 0.0;
    beta[7] = 0.0;
    const double doseThres(8.0);
    const double arrestTime(0.0);
    const int oxy(1);
    const double hypNecThres(0.0);
    double DO2(1.84);
    double VmaxO2(0.0152);
    const double KmO2(3.04);
    const double pO2NormVes(42.0);
    const double pO2TumVes(42.0);
    const double hypThres(5.0);


    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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
    treatment = new Treatment(x[0], x[0] * x[1], 24.0, 0);
    cout << treatment << endl;

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2, KmO2,
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

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This functions evaluates the reduced model of tumour growth and response to
 * radiotherapy to obtain its scalar outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void reducedModel(const double *x, double *y, const int nrow, const int ncol,
                  const int nlayer, const double cellSize,
                  const vector<bool> &inTum, const vector<bool> &inVes){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);

    const int edgeOrder(x[k]);
    k++;
    const double tumTime(x[k]);
    k++;
    vector<double> alpha(8), beta(8);
    alpha[1] = x[k];
    k++;
    beta[1]  = alpha[1] / x[k];
    k++;
    alpha[2] = x[k];
    k++;
    beta[2]  = alpha[2] / x[k];
    k++;
    alpha[3] = x[k];
    k++;
    beta[3]  = alpha[3] / x[k];
    k++;
    alpha[4] = x[k];
    k++;
    beta[4]  = alpha[4] / x[k];
    k++;
    alpha[5] = x[k];
    k++;
    beta[5]  = alpha[5] / x[k];
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

    const bool tumGrowth(true);
    const bool res(false);
    const double fibTime(1130.0);
    const bool ang(false);
    const double vascTumTime(2400.0);
    double Dvegf(2.0);
    double VmaxVegf(6.38e-3);
    const double KmVegf(2.5);
    const double vegfThres(15.0);
    const double hypVegf(20.0);
    alpha[0] = 0.0;
    beta[0]  = 0.0;
    alpha[6] = 0.0;
    beta[6] = 0.0;
    alpha[7] = 0.0;
    beta[7] = 0.0;
    const double doseThres(8.0);
    const double arrestTime(0.0);
    const double hypThres(5.0);
    const int oxy(1);
    const double pO2TumVes(42.0);

    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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
    treatment = new Treatment(2.0, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2, KmO2,
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

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This functions evaluates the reduced model of tumour growth and response to
 * radiotherapy to obtain its scalar outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void reducedModel(const double *x, double *y, const int nrow, const int ncol,
                  const int nlayer, const double cellSize,
                  const vector<bool> &inTum, const vector<bool> &inVes,
                  Treatment *const treatment){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);

    const int edgeOrder(x[k]);
    k++;
    const double tumTime(x[k]);
    k++;
    vector<double> alpha(8), beta(8);
    alpha[1] = x[k];
    k++;
    beta[1]  = alpha[1] / x[k];
    k++;
    alpha[2] = x[k];
    k++;
    beta[2]  = alpha[2] / x[k];
    k++;
    alpha[3] = x[k];
    k++;
    beta[3]  = alpha[3] / x[k];
    k++;
    alpha[4] = x[k];
    k++;
    beta[4]  = alpha[4] / x[k];
    k++;
    alpha[5] = x[k];
    k++;
    beta[5]  = alpha[5] / x[k];
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

    const bool tumGrowth(true);
    const bool res(false);
    const double fibTime(1130.0);
    const bool ang(false);
    const double vascTumTime(2400.0);
    double Dvegf(2.0);
    double VmaxVegf(6.38e-3);
    const double KmVegf(2.5);
    const double vegfThres(15.0);
    const double hypVegf(20.0);
    alpha[0] = 0.0;
    beta[0]  = 0.0;
    alpha[6] = 0.0;
    beta[6] = 0.0;
    alpha[7] = 0.0;
    beta[7] = 0.0;
    const double doseThres(8.0);
    const double arrestTime(0.0);
    const double hypThres(5.0);
    const int oxy(1);
    const double pO2TumVes(42.0);

    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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
    cout << treatment;

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2, KmO2,
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

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);

    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This functions evaluates the reduced model of tumour growth and response to
 * radiotherapy to obtain its scalar outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void reducedArtModel(const double *x, double *y){
    const double cellSize(20.0);
    const double tumArea(x[0]);
    const double tumDens(x[1]);
    const double vascDens(x[2]);
    int nrow, ncol, nlayer;
    vector<bool> inTum;
    vector<bool> inVes;

    createInFiles(cellSize, tumArea, tumDens, vascDens, nrow, ncol, nlayer,
                  inTum, inVes);

    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(3);
    const int edgeOrder(x[k]);
    k++;
    const double tumTime(x[k]);
    k++;
    vector<double> alpha(8), beta(8);
    alpha[1] = x[k];
    k++;
    beta[1]  = alpha[1] / x[k];
    k++;
    alpha[2] = x[k];
    k++;
    beta[2]  = alpha[2] / x[k];
    k++;
    alpha[3] = x[k];
    k++;
    beta[3]  = alpha[3] / x[k];
    k++;
    alpha[4] = x[k];
    k++;
    beta[4]  = alpha[4] / x[k];
    k++;
    alpha[5] = x[k];
    k++;
    beta[5]  = alpha[5] / x[k];
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

    const bool tumGrowth(true);
    const bool res(false);
    const double fibTime(1130.0);
    const bool ang(false);
    const double vascTumTime(2400.0);
    double Dvegf(2.0);
    double VmaxVegf(6.38e-3);
    const double KmVegf(2.5);
    const double vegfThres(15.0);
    const double hypVegf(20.0);
    alpha[0] = 0.0;
    beta[0]  = 0.0;
    alpha[6] = 0.0;
    beta[6] = 0.0;
    alpha[7] = 0.0;
    beta[7] = 0.0;
    const double doseThres(8.0);
    const double arrestTime(0.0);
    const double hypThres(5.0);
    const int oxy(1);
    const double pO2TumVes(42.0);

    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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
    treatment = new Treatment(2.0, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2, KmO2,
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

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This functions evaluates the reduced model of tumour growth and response to
 * radiotherapy to obtain both its scalar outputs and time dependent outputs.
 *
 * Inputs:
 *  - x: array containing the parameters of the model,
 *  - nFTumDens: name of the file where the tumour density values will be
 *  written,
 *  - nFTumVol: name of the file where the tumour volume values will be written,
 *  - nFVascDens: name of the file where the vascular density values will be
 *  written,
 *  - nFKilledCells: name of the file where the killed tumour cell percentages
 *  values will be written,
 *  - nFDeadDens: name of the file where the dead cell density values will be
 *  written,
 *  - nFCycle: name of the file where the cycle distributon values will be
 *  written,
 *  - nFHypDens: name of the file where the hypoxic density values will be
 *  written,
 *  - nFPO2Stat: name of the file where the pO2 statistic values will be
 *  written,
 *  - nFVegfStat: name of the file where the VEGF statistic values will be
 *  written.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void reducedArtModel(const double *x, double *y, const string nFTumDens,
                     const string nFTumVol, const string nFVascDens,
                     const string nFKilledCells, const string nFDeadDens,
                     const string nFCycle, const string nFHypDens,
                     const string nFPO2Stat, const string nFVegfStat){
    const double cellSize(20.0);
    const double tumArea(x[0]);
    const double tumDens(x[1]);
    const double vascDens(3.8);
    int nrow, ncol, nlayer;
    vector<bool> inTum;
    vector<bool> inVes;

    createInFiles(cellSize, tumArea, tumDens, vascDens, nrow, ncol, nlayer,
                  inTum, inVes);

    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(2);
    const int edgeOrder(x[k]);
    k++;
    const double tumTime(x[k]);
    k++;
    vector<double> alpha(8), beta(8);
    alpha[1] = x[k];
    k++;
    beta[1]  = alpha[1] / x[k];
    k++;
    alpha[2] = x[k];
    k++;
    beta[2]  = alpha[2] / x[k];
    k++;
    alpha[3] = x[k];
    k++;
    beta[3]  = alpha[3] / x[k];
    k++;
    alpha[4] = x[k];
    k++;
    beta[4]  = alpha[4] / x[k];
    k++;
    alpha[5] = x[k];
    k++;
    beta[5]  = alpha[5] / x[k];
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
    const bool tumGrowth(true);
    const bool res(false);
    const double fibTime(1130.0);
    const bool ang(false);
    const double vascTumTime(2400.0);
    double Dvegf(2.0);
    double VmaxVegf(6.38e-3);
    const double KmVegf(2.5);
    const double vegfThres(15.0);
    const double hypVegf(20.0);
    alpha[0] = 0.0;
    beta[0]  = 0.0;
    alpha[6] = 0.0;
    beta[6] = 0.0;
    alpha[7] = 0.0;
    beta[7] = 0.0;
    const double doseThres(8.0);
    const double arrestTime(0.0);
    const double hypThres(5.0);
    const int oxy(1);
    const double pO2TumVes(42.0);

    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
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
    treatment = new Treatment(2.0, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;
    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, edgeOrder,
                        tumGrowth, tumTime, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double oxySimTimeStep(10.0);
    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;
    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes, ang, Dvegf,
                           VmaxVegf, KmVegf, hypVegf, oxy, DO2, VmaxO2, KmO2,
                           pO2NormVes, pO2TumVes, hypThres);
    coupler = new Coupler(model1, model2);
    const double simTimeStep(6.0);
    const double sclFac(3.6e6 * simTimeStep / oxySimTimeStep);
    const double simTime(2160.0);

    RootSimulator *sim;
    sim = new RootSimulator(coupler, simTimeStep, oxySimTimeStep, sclFac);
    double currentTime(0.0);

    ofstream fTumDens(nFTumDens), fTumVol(nFTumVol);
    ofstream fKilledCells(nFKilledCells), fVascDens(nFVascDens);
    ofstream fDeadDens(nFDeadDens), fCycle(nFCycle), fHypDens(nFHypDens);
    ofstream fPO2Stat(nFPO2Stat), fVegfStat(nFVegfStat);

    sim->initSim();
    fTumDens     << currentTime << " " << model1->getOutD()[24] << endl;
    fTumVol      << currentTime << " " << model1->getOutD()[25] << endl;
    fVascDens    << currentTime << " " << model1->getOutD()[26] <<
                    " " << model1->getOutD()[27] <<
                    " " << model1->getOutD()[28] << endl;
    fKilledCells << currentTime << " " << model1->getOutD()[29] << endl;
    fHypDens     << currentTime << " " << model2->getOutD()[0] << endl;
    fDeadDens    << currentTime << " " << model1->getOutD()[30] << endl;
    fPO2Stat     << currentTime << " " << model2->getOutD()[1] << " " <<
                    model2->getOutD()[2] << endl;
    fVegfStat    << currentTime << " " << model2->getOutD()[3] << " " <<
                    model2->getOutD()[4] << endl;
    fCycle       << currentTime << " " << model1->getOutD()[31] <<
                    " " << model1->getOutD()[32] <<
                    " " << model1->getOutD()[33] <<
                    " " << model1->getOutD()[34] <<
                    " " << model1->getOutD()[35] << endl;

    int numIter(simTime / simTimeStep);
    for(int j(0); j < numIter; j++){
        currentTime += simTimeStep;
        sim->simulate(currentTime, simTimeStep);
        fTumDens     << currentTime << " " << model1->getOutD()[24] << endl;
        fTumVol      << currentTime << " " << model1->getOutD()[25] << endl;
        fVascDens    << currentTime << " " << model1->getOutD()[26] <<
                        " " << model1->getOutD()[27] <<
                        " " << model1->getOutD()[28] << endl;
        fKilledCells << currentTime << " " << model1->getOutD()[29] << endl;
        fHypDens     << currentTime << " " << model2->getOutD()[0]  << endl;
        fDeadDens    << currentTime << " " << model1->getOutD()[30] << endl;
        fPO2Stat     << currentTime << " " << model2->getOutD()[1] << " " <<
                        model2->getOutD()[2] << endl;
        fVegfStat    << currentTime << " " << model2->getOutD()[3] << " " <<
                        model2->getOutD()[4] << endl;
        fCycle       << currentTime << " " << model1->getOutD()[31] <<
                        " " << model1->getOutD()[32] <<
                        " " << model1->getOutD()[33] <<
                        " " << model1->getOutD()[34] <<
                        " " << model1->getOutD()[35] << endl;
    }

    fTumDens.close();
    fTumVol.close();
    fVascDens.close();
    fKilledCells.close();
    fDeadDens.close();
    fHypDens.close();
    fPO2Stat.close();
    fVegfStat.close();

    sim->stop();

    double eightWTumDens(model1->getOutD()[0]);
    double twelveWTumDens(model1->getOutD()[1]);
    double eightWTumVol(model1->getOutD()[2]);
    double twelveWTumVol(model1->getOutD()[3]);
    double eightWIntTumDens(model1->getOutD()[4]);
    double twelveWIntTumDens(model1->getOutD()[5]);
    double eightWIntTumVol(model1->getOutD()[6]);
    double twelveWIntTumVol(model1->getOutD()[7]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[11]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[12]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[21]);
    double recTumVol(model1->getOutD()[22]);
    double recTime(model1->getOutD()[23]);
    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "8wTumDens: "     << eightWTumDens  << endl;
    cout << "12wTumDens: "    << twelveWTumDens << endl;
    cout << "8wTumVol: "      << eightWTumVol  << " mm3" << endl;
    cout << "12wTumVol: "     << twelveWTumVol << " mm3" << endl;
    cout << "8wIntTumDens: "  << eightWIntTumDens << endl;
    cout << "12wIntTumDens: " << twelveWIntTumDens << endl;
    cout << "8wIntTumVol: "   << eightWIntTumVol << " mm3" << endl;
    cout << "12wIntTumVol: "  << twelveWIntTumVol << " mm3" << endl;
    cout << "killed50: "      << killed50 << endl;
    cout << "killed80: "      << killed80 << endl;
    cout << "killed90: "      << killed90 << endl;
    cout << "killed95: "      << killed95 << endl;
    cout << "timeTo95: "      << timeTo95 << " h" << endl;
    cout << "killed99: "      << killed99 << endl;
    cout << "timeTo99: "      << timeTo99 << " h" << endl;
    cout << "killed999: "     << killed999 << endl;
    cout << "rec: "           << rec << endl;
    cout << "recTumDens: "    << recTumDens << endl;
    cout << "recTumVol: "     << recTumVol << " mm3" << endl;
    cout << "recTime: "       << recTime << " h" << endl;

    y[0]  = eightWTumDens;
    y[1]  = twelveWTumDens;
    y[2]  = eightWTumVol;
    y[3]  = twelveWTumVol;
    y[4]  = eightWIntTumDens;
    y[5]  = twelveWIntTumDens;
    y[6]  = eightWIntTumVol;
    y[7]  = twelveWIntTumVol;
    y[8]  = killed50;
    y[9]  = killed80;
    y[10] = killed90;
    y[11] = killed95;
    y[12] = timeTo95;
    y[13] = killed99;
    y[14] = timeTo99;
    y[15] = killed999;
    y[16] = rec;
    y[17] = recTumDens;
    y[18] = recTumVol;
    y[19] = recTime;
}


/*------------------------------------------------------------------------------
 * This function evaluates the toy model.
 *
 * Inputs:
 *  - x: array containing the parameters of the model.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void toyModel(double *x, double *y){
    y[0] = x[0] + 2* x[1] + x[2] * x[2] + x[3] * x[4];
}
