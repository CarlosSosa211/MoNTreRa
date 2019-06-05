#include "evalModel.hpp"

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
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - inPO2: vector containing the initial pO2 values.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void model(const double *x, double *y, const int nrow, const int ncol,
           const int nlayer, const double cellSize, const vector<bool> &inTum,
           const vector<bool> &inVes, const vector<double> &inPO2){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);
    const int tumGrowth(x[k] > 0.5);
    k++;
    const double tumTime(x[k]);
    k++;
    const int edgeOrder(x[k]);
    k++;
    const int res(x[k] > 0.5);
    k++;
    const double fibTime(x[k]);
    k++;
    const int ang(x[k] > 0.5);
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

    double endTreatTumDens(model1->getOutD()[0]);
    double threeMonTumDens(model1->getOutD()[1]);
    double tumVol(model1->getOutD()[19]);
    double intTumDens(model1->getOutD()[2]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[6]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[7]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[16]);
    double recTime(model1->getOutD()[17]);

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
    const int ang(x[k] > 0.5);
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

    std::ofstream fPO2(nFPO2.c_str());

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
           const vector<bool> &inVes, const vector<double> &inPO2,
           const string nFTumDens, const string nFTumVol,
           const string nFVascDens, const string nFKilledCells,
           const string nFDeadDens, const string nFCycle,
           const string nFHypDens, const string nFPO2Stat,
           const string nFVegfStat){
    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};

    int k(0);
    const int tumGrowth(x[k]);
    k++;
    const double tumTime(x[k]);
    k++;
    const int edgeOrder(x[k]);
    k++;
    const int res(x[k]);
    k++;
    const double fibTime(x[k]);
    k++;
    const int ang(x[k]);
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
    AbsOxyTissue *model2;

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


    switch(oxy){
    case 0:{
        const std::vector<double> inPO2(nrow * ncol * nlayer, 0.0);
        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, inPO2,
                                    oxy, hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = 1.0;
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
        break;
    }

    case 4:{
        const int nComp(nrow * ncol * nlayer);
        std::vector<double> inPO2(nComp, 0.0);

        model2 = new ConstOxyTissue(nrow, ncol, nlayer, inVes, inPO2, oxy,
                                    hypThres);
        coupler = new Coupler(model1, model2);
        sclFac = 1.0;
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
    ofstream fPO2Stat(nFPO2Stat), fVEGFStat(nFVegfStat);

    sim->initSim();

    fTumDens     << currentTime << " " << model1->getOutD()[18] << std::endl;
    fTumVol      << currentTime << " " << model1->getOutD()[19] << std::endl;
    fVascDens    << currentTime << " " << model1->getOutD()[20] <<
                    " " << model1->getOutD()[21] <<
                    " " << model1->getOutD()[22] << std::endl;
    fKilledCells << currentTime << " " << model1->getOutD()[23] << std::endl;
    fHypDens     << currentTime << " " << model2->getOutD()[0] << std::endl;
    fDeadDens    << currentTime << " " << model1->getOutD()[24] << std::endl;
    fPO2Stat     << currentTime << " " << model2->getOutD()[1] << " " <<
                    model2->getOutD()[2] << std::endl;
    fVEGFStat    << currentTime << " " << model2->getOutD()[3] << " " <<
                    model2->getOutD()[4] << std::endl;
    fCycle       << currentTime << " " << model1->getOutD()[25] <<
                    " " << model1->getOutD()[26] <<
                    " " << model1->getOutD()[27] <<
                    " " << model1->getOutD()[28] <<
                    " " << model1->getOutD()[29] << std::endl;

    int numIter(simTime / simTimeStep);

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

    double endTreatTumDens(model1->getOutD()[0]);
    double threeMonTumDens(model1->getOutD()[1]);
    double tumVol(model1->getOutD()[19]);
    double intTumDens(model1->getOutD()[2]);
    double killed50(model1->getOutB()[0]);
    double killed80(model1->getOutB()[1]);
    double killed90(model1->getOutB()[2]);
    double killed95(model1->getOutB()[3]);
    double timeTo95(model1->getOutD()[6]);
    double killed99(model1->getOutB()[4]);
    double timeTo99(model1->getOutD()[7]);
    double killed999(model1->getOutB()[5]);
    double rec(model1->getOutB()[6]);
    double recTumDens(model1->getOutD()[16]);
    double recTime(model1->getOutD()[17]);

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
         const string nFInTissueDim, const string nFInVes,
         const string nFInPO2){
    const int K(12), nOut(1);
    bool art(0);
    int nrow, ncol, nlayer;
    double cellSize, vascDens, sigmaVasc;
    vector<bool> inVes;
    vector<double> inPO2;

    readInFilesOxy(nFInTissueOxy, art, nrow, ncol, nlayer, cellSize, vascDens,
                   sigmaVasc);

    if(!art){
        readInFiles(nFInTissueDim, nFInVes, nFInPO2, nrow, ncol, nlayer,
                    cellSize, inVes, inPO2);
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
 * This function evaluates the toy model.
 *
 * Inputs:
 *  - x: array containing the parameters of the model.
 *
 * Outputs:
 *  - y: array containing the scalar outputs of the model.
------------------------------------------------------------------------------*/

void toyModel(double *x, double *y){
    y[0] = x[0] + 2.0 * x[1] + x[2] * x[2] + x[3] * x[4];
}
