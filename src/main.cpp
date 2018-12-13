#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "oxyTissue.hpp"
#include "rootSimulator.hpp"
#include "sensAn.hpp"
#include "simpcell.hpp"
#include "tissue.hpp"
#include "treatment.hpp"

using namespace std;


int main(){
    //const int N(1e5);
    const int kp(5), L(2), p(20), N(100);
    //const int nMethod(1), nModel(0);
    //string nFRefParInt("../InputFiles/refParIntOneAlphaBeta.dat");
    //string nFRefParInt("../InputFiles/refParIntRTAll.dat");
    //string nFRefParInt("../InputFiles/refParIntToy.dat");
    //string nRefParMean("../InputFiles/refParMeanRT.dat");
    //string nFMostRelPar("../InputFiles/mostRelParAng.dat");
    //string nFLeastRelPar("../InputFiles/leastRelParAng.dat");
    //string nFVarPar("../InputFiles/varParAng.dat");
    string nFMostRelPar("../InputFiles/mostRelParAlphaBeta.dat");
    string nFLeastRelPar("../InputFiles/leastRelParAlphaBeta.dat");
    string nFVarPar("../InputFiles/varParAlphaBeta.dat");
    //srand(time(NULL));
    //var1ParRange(kp, L, nFRefParInt);
    varErr(nFVarPar, nFMostRelPar, nFLeastRelPar);
    //vector<string> nFPar;
    //nFPar.push_back("../InputFiles/parDiffAlphaBeta.dat");
    //nFPar.push_back("../InputFiles/parOneAlphaBeta.dat");
    //varParFromFiles(nFPar);
    //evalR(nMethod, nModel);
    //morrisRT(N, p, nFRefParInt);
    //morrisToy(N, p, nFRefParInt);
    //morrisVarRangeRT(kp, L, N, p, nFRefParInt);
    //morrisVarRangeToy(kp, L, N, p, nFRefParInt);
    //sobolRT(N, nFRefParInt);
    //sobolToy(N, nFRefParInt);
    //sobolFromFiles(2);

    return 0;
}


int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double sigmaTum,
                  const double vascDens, const double sigmaVasc,
                  vector<bool> &inTum, vector<bool> &inVes){
    int ivd, ivd2;
    int mindim, mindim2, sqrtmin, tumToDist, vesToDist;
    int nrowNcol, nrowNcolNlayer;
    vector<int> div;
    vector<double> diff;

    nrowNcol = nrow * ncol;
    nrowNcolNlayer = nrowNcol * nlayer;
    tumToDist = tumDens * nrowNcolNlayer;

    if(vascDens){
        mindim = min(nrow, ncol);
        mindim2 = mindim * mindim;
        sqrtmin = sqrt(mindim);

        for(int l(1); l < sqrtmin; l+=2){
            if(!(nrow % l) && !(ncol % l)){
                div.push_back(l);
                diff.push_back(fabs(1.0 / (l * l) - vascDens));
                div.push_back(mindim / l);
                diff.push_back(fabs(double(l * l) / double(mindim2) - vascDens));
            }
        }

        ivd = div.at(min_element(diff.begin(), diff.end()) - diff.begin());
        ivd2 = ivd * ivd;
        vesToDist = min(nrowNcol / ivd2, nrowNcolNlayer - tumToDist);
    }

    else{
        vesToDist = 0;
    }

    int halfIvd, halfNrow, halfNcol, halfNlayer;
    int imHalfNrow2, lmHalfNlayer2;
    int lnrowNcol;
    int irr2, ishift;
    double R, RR;
    vector<simpCell> map(nrowNcolNlayer);

    halfIvd = 0.5 * ivd;
    halfNrow = 0.5 * nrow;
    halfNcol = 0.5 * ncol;
    halfNlayer = 0.5 * nlayer;
    R = sqrt(halfNrow * halfNrow + halfNcol * halfNcol + halfNlayer * halfNlayer);
    RR = 1.001 * sqrt(2.0 * halfIvd * halfIvd);

    int k(0);

    for(int l(0); l < nlayer; l++){
        lmHalfNlayer2 = (l - halfNlayer) * (l - halfNlayer);
        lnrowNcol = l * nrowNcol;
        for(int i(0); i < nrow; i++){
            imHalfNrow2 = (i - halfNrow) * (i - halfNrow);
            irr2 = (i % ivd - halfIvd) * (i % ivd - halfIvd);
            ishift = i / ivd * ncol / ivd;
            for(int j(0); j < ncol; j++){
                map.at(k).r = sqrt(lmHalfNlayer2 + imHalfNrow2 +
                                   (j - halfNcol) * (j - halfNcol)) / R;
                map.at(k).rr = lnrowNcol + ishift + j / ivd +
                        sqrt(irr2 + (j % ivd - halfIvd) * (j % ivd - halfIvd)) / RR;
                map.at(k).k = k;
                map.at(k).tum = 0;
                map.at(k).ves = 0;
                k++;
            }
        }
    }
    k = 0;

    sort(map.begin(), map.end(), compR);

    int m;
    double n;
    default_random_engine gen;
    normal_distribution<double> distTum(0, sigmaTum);
    bool tooBigTumDens(tumDens > 0.9);
    bool tooSmallSigmaTum(sigmaTum < 0.15);

    if(!tooBigTumDens && !tooSmallSigmaTum){
        while(tumToDist > 0){
            n = distTum(gen);
            if(n >= 0.0 && n < 1.0){
                m = n * nrowNcolNlayer;
                if(!map.at(m).tum){
                    map.at(m).tum = 1;
                    tumToDist--;
                }
            }
        }
    }

    else if(tooBigTumDens){
        for(int k(0); k < map.size(); k++){
            map.at(k).tum = 1;
        }
        reverse(map.begin(), map.end());
        int tumToRemove(nrowNcolNlayer - tumToDist);
        while(tumToRemove > 0){
            n = distTum(gen);
            if(n >= 0.0 && n < 1.0){
                m = n * nrowNcolNlayer;
                if(map.at(m).tum){
                    map.at(m).tum = 0;
                    tumToRemove--;
                }
            }
        }
    }

    else{
        bool cond;
        while(tumToDist > 0){
            n = distTum(gen);
            if(n >= 0.0 && n < 1.0){
                m = n * nrowNcolNlayer;
                cond = true;
                while(cond){
                    if(!map.at(m).tum){
                        map.at(m).tum = 1;
                        tumToDist--;
                        cond = false;
                    }
                    else{
                        m++;
                    }
                    if(m == nrowNcolNlayer){
                        m = 0;
                    }
                }
            }
        }
    }

    sort(map.begin(), map.end(), compRR);

    int nivd2;
    bool availSpace, firstRound(true);
    normal_distribution<double> distVes(0, sigmaVasc);

    k = 0;
    while(vesToDist > 0){
        availSpace = false;
        for(int j(0); j < ivd2; j++){
            if(!map.at(k + j).ves && !map.at(k + j).tum){
                availSpace = true;
                break;
            }
        }

        if(availSpace){
            n = distVes(gen);
            while(n < 0.0 || n >= 1.0){
                n = distVes(gen);
            }
            nivd2 = n * ivd2;
            m = k + nivd2;
            while(map.at(m).ves || map.at(m).tum){
                nivd2++;
                if(nivd2 >= ivd2){
                    nivd2 = 0;
                }
                m = k + nivd2;
            }
            map.at(m).ves = 1;
            vesToDist--;
            availSpace = false;
        }

        if(firstRound){
            k += ivd2;
        }
        else{
            k = rand() % (nrowNcol / ivd2) * ivd2;
        }
        if(k >= nrowNcol - ivd2){
            firstRound = false;
        }
    }

    sort(map.begin(), map.end(), compK);
    for(int k(0); k < nrowNcol; k++){
        if(map.at(k).ves == 1){
            for(int l(1); l < nlayer; l++){
                map.at(l * nrowNcol + k).ves = 1;
                map.at(l * nrowNcol + k).tum = 0;
            }
        }
    }

    for(k = 0; k < nrowNcolNlayer; k++){
        inTum.push_back(map.at(k).tum);
        inVes.push_back(map.at(k).ves);
    }

    return 0;
}


void model(const double *x, double *y){
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

    /*int nrow, ncol, nlayer;
    double cellSize;

    std::ifstream fInTissueDim("../InputFiles/tissueDim.dat");

    fInTissueDim >> nrow >> ncol >> nlayer;
    fInTissueDim >> cellSize;
    fInTissueDim.close();

    vector<bool> inTum;
    std::ifstream fInTum("../InputFiles/inTum.dat");
    bool temp;

    fInTum >> temp;
    while(!fInTum.eof()){
        inTum.push_back(temp);
        fInTum >> temp;
    }
    fInTum.close();

    vector<bool> inVes;
    std::ifstream fInVes("../InputFiles/inVes.dat");

    fInVes >> temp;
    while(!fInVes.eof()){
        inVes.push_back(temp);
        fInVes >> temp;
    }
    fInVes.close();*/

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

    int nrow, ncol, nlayer;
    double cellSize;
    std::ifstream fInTissueDim("../InputFiles/tissueDim.dat");

    fInTissueDim >> nrow >> ncol >> nlayer;
    fInTissueDim >> cellSize;
    fInTissueDim.close();

    vector<bool> inTum;
    std::ifstream fInTum("../InputFiles/inTum.dat");
    bool temp;

    fInTum >> temp;
    while(!fInTum.eof()){
        inTum.push_back(temp);
        fInTum >> temp;
    }
    fInTum.close();

    vector<bool> inVes;
    std::ifstream fInVes("../InputFiles/inVes.dat");

    fInVes >> temp;
    while(!fInVes.eof()){
        inVes.push_back(temp);
        fInVes >> temp;
    }
    fInVes.close();

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
    const double vegfThres(x[k]);
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

    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "res: "         << res         << endl;
    cout << "fibTime: "     << fibTime     << " h" << endl;
    cout << "ang: "         << ang         << endl;
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
    cout << "hypVegf: "     << hypVegf     << " mol/um^3" << endl;

    Treatment *treatment;
    treatment = new Treatment(dose, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes,
                        tumGrowth, tumTime, edgeOrder, cycDur, cycDistrib,
                        res, fibTime, ang, vascTumTime, vegfThres, alpha,
                        beta, doseThres, arrestTime, treatment, hypNecThres);

    double simTimeStep, oxySimTimeStep, sclFac, simTime;

    simTimeStep    = 6.0; //h;
    oxySimTimeStep = 10.0; //ms;
    sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;

    DO2      *= oxySimTimeStep;
    Vmax     *= oxySimTimeStep;
    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes,
                           Dvegf, DO2, Vmax, Km, pO2NormVes, pO2TumVes,
                           hypThres, ang, VmaxVegf, KmVegf, hypVegf);

    coupler = new Coupler(model1, model2);

    RootSimulator *sim;

    sim = new RootSimulator(coupler, simTimeStep,
                            oxySimTimeStep, sclFac);
    simTime = treatment->getDuration() + 720.0; // + 3 months

    sim->initSim();
    sim->simulate(simTimeStep, simTime);
    sim->stop();

    double endTreatTumDens(model1->getOut()->at(24));
    double threeMonTumDens(model1->getOut()->at(25));
    double recTumDens(model1->getOut()->at(26));
    double tumVol(model1->getOut()->at(23));
    double intTumDens(model1->getOut()->at(22));
    double timeTo95(model1->getOut()->at(12));
    double timeTo99(model1->getOut()->at(13));
    double recTime(model1->getOut()->at(27));

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "endTreatTumDens: " << endTreatTumDens << endl;
    cout << "3MonTumDens: "     << threeMonTumDens << endl;
    cout << "recTumDens: "      << recTumDens << endl;
    cout << "tumVol: "          << tumVol << " mm3" << endl;
    cout << "intTumDens: "      << intTumDens << endl;
    cout << "timeTo95: "        << timeTo95 << " h" << endl;
    cout << "timeTo99: "        << timeTo99 << " h" << endl;
    cout << "recTime: "         << recTime << " h" << endl;

    y[0] = endTreatTumDens;
    y[1] = threeMonTumDens;
    y[2] = recTumDens;
    y[3] = tumVol;
    y[4] = intTumDens;
    y[5] = timeTo95;
    y[6] = timeTo99;
    y[7] = recTime;
}


void model(const double *x, double *y, const std::string nFTumDens,
           const std::string nFTumVol, const std::string nFVascDens,
           const std::string nFKilledCells, const std::string nFCycle,
           const std::string nFHypDens, const std::string nFPO2Stat,
           const std::string nFVegfStat){
    int nrow, ncol, nlayer;
    double cellSize;
    std::ifstream fInTissueDim("../InputFiles/tissueDim.dat");

    fInTissueDim >> nrow >> ncol >> nlayer;
    fInTissueDim >> cellSize;
    fInTissueDim.close();

    vector<bool> inTum;
    std::ifstream fInTum("../InputFiles/inTum.dat");
    bool temp;

    fInTum >> temp;
    while(!fInTum.eof()){
        inTum.push_back(temp);
        fInTum >> temp;
    }
    fInTum.close();

    vector<bool> inVes;
    std::ifstream fInVes("../InputFiles/inVes.dat");

    fInVes >> temp;
    while(!fInVes.eof()){
        inVes.push_back(temp);
        fInVes >> temp;
    }
    fInVes.close();

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
    const double vegfThres(x[k]);
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

    cout << "tumGrowth: "   << tumGrowth   << endl;
    cout << "tumTime: "     << tumTime     << " h" << endl;
    cout << "edgeOrder: "   << edgeOrder   << endl;
    cout << "res: "         << res         << endl;
    cout << "fibTime: "     << fibTime     << " h" << endl;
    cout << "ang: "         << ang         << endl;
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
    cout << "hypVegf: "     << hypVegf     << " mol/um^3" << endl;

    Treatment *treatment;
    treatment = new Treatment(dose, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes,
                        tumGrowth, tumTime, edgeOrder, cycDur, cycDistrib,
                        res, fibTime, ang, vascTumTime, vegfThres, alpha,
                        beta, doseThres, arrestTime, treatment, hypNecThres);

    const double simTimeStep(6.0), oxySimTimeStep(10.0);

    DO2      *= oxySimTimeStep;
    Vmax     *= oxySimTimeStep;
    Dvegf    *= oxySimTimeStep;
    VmaxVegf *= oxySimTimeStep;

    model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes,
                           Dvegf, DO2, Vmax, Km, pO2NormVes, pO2TumVes,
                           hypThres, ang, VmaxVegf, KmVegf, hypVegf);

    coupler = new Coupler(model1, model2);

    const double sclFac(3.6e6 * simTimeStep / oxySimTimeStep);
    const double simTime(treatment->getDuration() + 720.0);
    RootSimulator *sim;

    sim = new RootSimulator(coupler, simTimeStep,
                            oxySimTimeStep, sclFac);

    int numIter(simTime / simTimeStep);
    double currentTime(0.0);
    std::ofstream fTumDens(nFTumDens), fTumVol(nFTumVol), fKilledCells(nFKilledCells);
    std::ofstream fVascDens(nFVascDens), fCycle(nFCycle), fHypDens(nFHypDens);
    std::ofstream fPO2Stat(nFPO2Stat), fVEGFStat(nFVegfStat);

    sim->initSim();

    fTumDens     << currentTime << " " <<
                    model1->getOut()->at(0) << std::endl;
    fTumVol      << currentTime << " " <<
                    model1->getOut()->at(23) << std::endl;
    fVascDens    << currentTime << " " <<
                    model1->getOut()->at(6) << " " <<
                    model1->getOut()->at(7) << " " <<
                    model1->getOut()->at(8) << std::endl;
    fKilledCells << currentTime << " " <<
                    model1->getOut()->at(21) << std::endl;
    fHypDens     << currentTime << " " <<
                    coupler->getModel2()->getOut()->at(0) << std::endl;
    fPO2Stat     << currentTime << " " <<
                    coupler->getModel2()->getOut()->at(1) << " " <<
                    coupler->getModel2()->getOut()->at(2) << std::endl;
    fVEGFStat    << currentTime << " " <<
                    coupler->getModel2()->getOut()->at(3) << " " <<
                    coupler->getModel2()->getOut()->at(4) << std::endl;
    fCycle       << currentTime << " " <<
                    model1->getOut()->at(1) << " " <<
                    model1->getOut()->at(2) << " " <<
                    model1->getOut()->at(3) << " " <<
                    model1->getOut()->at(4) << " " <<
                    model1->getOut()->at(5) << std::endl;

    for(int j(0); j < numIter; j++){
        currentTime += simTimeStep;
        sim->simulate(currentTime, simTimeStep);

        fTumDens     << currentTime << " " <<
                        model1->getOut()->at(0) << std::endl;
        fTumVol      << currentTime << " " <<
                        model1->getOut()->at(23) << std::endl;
        fVascDens    << currentTime << " " <<
                        model1->getOut()->at(6) << " " <<
                        model1->getOut()->at(7) << " " <<
                        model1->getOut()->at(8) << std::endl;
        fKilledCells << currentTime << " " <<
                        model1->getOut()->at(21) << std::endl;
        fHypDens     << currentTime << " " <<
                        coupler->getModel2()->getOut()->at(0) << std::endl;
        fPO2Stat     << currentTime << " " <<
                        coupler->getModel2()->getOut()->at(1) << " " <<
                        coupler->getModel2()->getOut()->at(2) << std::endl;
        fVEGFStat    << currentTime << " " <<
                        coupler->getModel2()->getOut()->at(3) << " " <<
                        coupler->getModel2()->getOut()->at(4) << std::endl;
        fCycle       << currentTime << " " <<
                        model1->getOut()->at(1) << " " <<
                        model1->getOut()->at(2) << " " <<
                        model1->getOut()->at(3) << " " <<
                        model1->getOut()->at(4) << " " <<
                        model1->getOut()->at(5) << std::endl;
    }

    fTumDens.close();
    fTumVol.close();
    fVascDens.close();
    fKilledCells.close();
    fHypDens.close();
    fPO2Stat.close();
    fVEGFStat.close();

    sim->stop();

    double endTreatTumDens(model1->getOut()->at(24));
    double threeMonTumDens(model1->getOut()->at(25));
    double recTumDens(model1->getOut()->at(26));
    double tumVol(model1->getOut()->at(23));
    double intTumDens(model1->getOut()->at(22));
    double timeTo95(model1->getOut()->at(12));
    double timeTo99(model1->getOut()->at(13));
    double recTime(model1->getOut()->at(27));

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "endTreatTumDens: " << endTreatTumDens << endl;
    cout << "3MonTumDens: "     << threeMonTumDens << endl;
    cout << "recTumDens: "      << recTumDens      << endl;
    cout << "tumVol: "          << tumVol          << " mm3" << endl;
    cout << "intTumDens: "      << intTumDens      << endl;
    cout << "timeTo95: "        << timeTo95        << " h" << endl;
    cout << "timeTo99: "        << timeTo99        << " h" << endl;
    cout << "recTime: "         << recTime         << " h" << endl;

    y[0] = endTreatTumDens;
    y[1] = threeMonTumDens;
    y[2] = recTumDens;
    y[3] = tumVol;
    y[4] = intTumDens;
    y[5] = timeTo95;
    y[6] = timeTo99;
    y[7] = recTime;
}


void toyModel(double *x, double *y){
    y[0] = x[0] + 2.0 * x[1] + x[2] * x[2] + x[3] * x[4];
}
