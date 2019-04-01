#include "tcp.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * This functions evaluates the model of tumour growth and response to
 * radiotherap considering a given treatment. Simulation is stopped when the
 * tumour is controlled.
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
 *  - treatment: pointer to the treatment to be considered
 *
 * Outputs:
 *  - y: array containing the scalar outputs "controlled" and "doseToControl" of
 *  the model.
------------------------------------------------------------------------------*/

void modelTCP(const double *x, double *y, const int nrow, const int ncol,
              const int nlayer, const double cellSize,
              const vector<bool> &inTum, const vector<bool> & inVes,
              Treatment *const treatment){
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

    cout << "tumGrowth: " << tumGrowth << endl;
    cout << "tumTime: "   << tumTime   << " h" << endl;
    cout << "edgeOrder: " << edgeOrder << endl;
    cout << "res: "       << res       << endl;
    if(res){
        cout << "fibTime: " << fibTime << " h" << endl;
    }
    cout << "ang: " << ang << endl;
    if(ang){
        cout << "vascTumTime: " << vascTumTime << " h" << endl;
        cout << "Dvegf: "       << Dvegf       << " um^2/ms" << endl;
        cout << "VmaxVegf: "    << VmaxVegf    << " mol/um^3ms" << endl;
        cout << "KmVegf: "      << KmVegf      << " mol/um^3" << endl;
        cout << "vegfThres: "   << vegfThres   << " mol/um^3" << endl;
        cout << "hypVegf: "     << hypVegf     << " mol/um^3" << endl;
    }
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
    cout << "doseThres: "  << doseThres  << " Gy" << endl;
    cout << "arrestTime: " << arrestTime << " h" << endl;
    cout << "oxy: "        << oxy        << endl;
    if(oxy){
        cout << "hypNecThres: " << hypNecThres << " mmHg" << endl;
        cout << "DO2: "         << DO2         << " um^2/ms" << endl;
        cout << "VmaxO2: "      << VmaxO2      << " mmHg/ms" << endl;
        cout << "KmO2: "        << KmO2        << " mmHg" << endl;
        cout << "pO2NormVes: "  << pO2NormVes  << " mmHg" << endl;
        cout << "pO2TumVes: "   << pO2TumVes   << " mmHg" << endl;
        cout << "hypThres: "    << hypThres    << " mmHg" << endl;
    }
    cout << treatment;

    double sclFac;
    Coupler *coupler;
    Tissue *model1;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes, tumGrowth,
                        tumTime, edgeOrder, cycDur, cycDistrib, res, fibTime,
                        ang, vascTumTime, vegfThres, alpha, beta, treatment,
                        doseThres, arrestTime, oxy, hypNecThres);

    const double simTimeStep(6.0);
    const double oxySimTimeStep(10.0);

    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;
    DO2      *= oxySimTimeStep;
    VmaxO2   *= oxySimTimeStep;

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

    const double simTime(treatment->getDuration());
    double currentTime(0.0);
    RootSimulator *sim;

    sim = new RootSimulator(coupler, simTimeStep, oxySimTimeStep, sclFac);

    sim->initSim();

    double controlled(0.0);

    const int numIter(simTime / simTimeStep);
    int j(0);

    while(j < numIter && !controlled){
        currentTime += simTimeStep;
        sim->simulate(currentTime, simTimeStep);
        controlled = model1->getOut()->at(35);
        j++;
    }

    sim->stop();

    const double doseToControl(model1->getOut()->at(36));

    delete model1;
    delete coupler->getModel2();
    delete coupler;
    delete sim;

    cout << "controlled: "    << controlled    << endl;
    cout << "doseToControl: " << doseToControl << " Gy" << endl;

    y[0] = controlled;
    y[1] = doseToControl;
}


/*------------------------------------------------------------------------------
 * This functions performs the simulations needed to obtain tcp curves.
 *
 * Inputs:
 *  - N: number of repetitions,
 *  - nFInTissueTCP: name of the file containing the architecture of the tissue
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
 *  - treatment: pointer to the treatment to be considered
 *
 * Outputs:
 *  - y: array containing the scalar outputs "controlled" and "doseToControl" of
 *  the model.
------------------------------------------------------------------------------*/

void tcp(const int N, const string nFInTissueTCP, const string nFParTCP,
         const vector<string> nFTreatmentTCP, const string nFInTissueDim,
         const string nFInTum, const string nFInVes){
    const int K(38), nOut(2);
    bool art(0);
    int nrow, ncol, nlayer;
    double cellSize, tumDens, sigmaTum, vascDens, sigmaVasc;
    vector<bool> inTum, inVes;
    vector<Treatment> treatment;

    readInFilesTCP(nFInTissueTCP, nFTreatmentTCP, art, nrow, ncol, nlayer,
                   cellSize, tumDens, sigmaTum, vascDens, sigmaVasc, treatment);
    if(!art){
        readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                    cellSize, inTum, inVes);
    }

    double x[K], y[nOut];
    ifstream fParTCP(nFParTCP.c_str());

    for(int k(0); k < K; k++){
        fParTCP >> x[k];
    }
    fParTCP.close();

    const int nEvTot(treatment.size() * N);
    int nEv(0);
    for(int i(0); i < treatment.size(); i++){
        ofstream fControlled("../OutputFiles/controlled_" + to_string(i) +
                             ".res");
        ofstream fDoseToControl("../OutputFiles/doseToControl_" + to_string(i) +
                                ".res");
        for(int j(0); j < N; j++){
            if(art){
                createInFiles(nrow, ncol, nlayer, tumDens, sigmaTum,
                              vascDens, sigmaVasc, inTum, inVes);
            }
            modelTCP(x, y, nrow, ncol, nlayer, cellSize, inTum, inVes,
                     &(treatment[i]));
            nEv++;
            cout << nEv << " out of " << nEvTot <<
                    " evaluations of the model" << endl;
            cout << "---------------------------------------------" << endl;

            fControlled    << y[0] << endl;
            fDoseToControl << y[1] << endl;
        }
        fControlled.close();
        fDoseToControl.close();
    }
}
