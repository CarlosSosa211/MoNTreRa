#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <time.h>
#include <vector>

#include "oxyTissue.hpp"
#include "rootSimulator.hpp"
#include "simpcell.hpp"
#include "tissue.hpp"
#include "treatment.hpp"

using namespace std;

int createInFiles(int nrow, int ncol, int nlayer,
                  double tumDens, double sigmaTum,
                  double vascDens, double sigmaVasc,
                  vector<bool> &inTum, vector<bool> &inVes);
void evalModelMorrisR();
void evalModelSobolR();
vector<double> model(vector<double> x);
void morris(int K, int p, int N);
void morrisk(int kp, int K, int p, int N);
void morrisR();
void sobol(int K, int N);
void sobolFromFiles(int K);
double toyModel(vector<double> x);

int main(){
    //srand(time(NULL));
    //int K(33), N(1e4);
    //sobol(K, N);

    //vector<double> x = {0.8, 0.8};
    //model(x);

    //sobolFromFiles(2);

     int K(34), p(20), N(100);
	 morris(K, p, N);
    //morrisk(0, K, p, N);
    //evalModelMorrisR();
    //evalModelSobolR();
    return 0;
}


void evalModelMorrisR(){
    ifstream fDim("../InputFiles/morrisDim.dat");
    int K, N;
    fDim >> K >> N;
    fDim.close();

    vector<double> x(K), y(4);
    ifstream fX("../InputFiles/X.dat");
    ofstream fYTumDens("../OutputFiles/YTumDens.res");
    ofstream fYIntTumDens("../OutputFiles/YIntTumDens.res");
    ofstream fYTimeTo95("../OutputFiles/YTimeTo95.res");
    ofstream fYTimeTo99("../OutputFiles/YTimeTo99.res");
    int nEv((K + 1) * N);
    for(int i(0); i < nEv; i++){
        for(int k(0); k < K; k++){
            fX >> x[k];
        }
	y = model(x);
    fYTumDens << y[0] << endl;
	fYIntTumDens << y[1] << endl;
	fYTimeTo95 << y[2] << endl;
	fYTimeTo99 << y[3] << endl;

        cout << i + 1 << " out of " << nEv << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
    }
    fX.close();
    fYTumDens.close();
    fYIntTumDens.close();
    fYTimeTo95.close();
    fYTimeTo99.close();
}


void evalModelSobolR(){
    ifstream fDim("../InputFiles/sobolDim.dat");
    int K, N;
    fDim >> K >> N;
    fDim.close();

    vector<double> x(K);
    ifstream fX("../InputFiles/X.dat");
    ofstream fY("../OutputFiles/Y.res");
    int nEv((K + 2) * N);
    for(int i(0); i < nEv; i++){
        for(int k(0); k < K; k++){
            fX >> x[k];
        }
        fY << toyModel(x) << endl;
        //cout << i + 1 << " out of " << nEv << " evaluations of the model" << endl;
        //cout << "---------------------------------------------" << endl;
    }
    fX.close();
    fY.close();
}


int createInFiles(int nrow, int ncol, int nlayer,
                  double tumDens, double sigmaTum,
                  double vascDens, double sigmaVasc,
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


vector<double> model(vector<double> x){
    int k(0);
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

    int nrow, ncol, nlayer;

    std::ifstream fInTissueDim("../InputFiles/tissueDim.dat");

    fInTissueDim >> nrow >> ncol >> nlayer;
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

    /*createInFiles(nrow, ncol, nlayer, tumDens, sigmaTum,
                  vascDens, sigmaVasc, inTum, inVes);*/

    vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    double tumGrowth(1.0);
    int edgeOrder(1);
    vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};
    double res(1.0);
    double ang(1.0);


    double tumTime(x[k]);
    k++;
    double fibTime(x[k]);
    k++;
    double vascTumTime(x[k]);
    k++;
    double vegfThres(x[k]);
    k++;
    vector<double> alpha(8), beta(8);
    for(int i(0); i < 8; i++){
        alpha[i] = x[k];
        k++;
        beta[i]  = alpha[i] / x[k];
        k++;
    }
    double doseThres(x[k]);
    k++;
    double arrestTime(x[k]);
    k++;
    double hypNecThres(x[k]);
    k++;
    double dose(x[k]);
    k++;
    double DO2(x[k]);
    k++;
    double Vmax(x[k]);
    k++;
    double Km(x[k]);
    k++;
    double Dvegf(x[k]);
    k++;
    double VmaxVegf(x[k]);
    k++;
    double KmVegf(x[k]);
    k++;
    double pO2NormVes(x[k]);
    k++;
    double pO2TumVes(x[k]);
    k++;
    double hypThres(x[k]);
    k++;
    double hypVegf(x[k]);
    k++;

    cout << "tumTime: " << tumTime << " h" <<endl;
    cout << "fibTime: " << fibTime << " h" << endl;
    cout << "vascTumTime: " << vascTumTime << " h" << endl;
    cout << "vegfThres: " << vegfThres << " mol/um^3" << endl;
    cout << "alpha : ";
    for(int i(0); i < 8; i++){
        cout << alpha[i] << " ";
    }
    cout << "Gy^-1" << endl;
    cout << "beta : ";
    for(int i(0); i < 8; i++){
        cout << beta[i] << " ";
    }
    cout << "Gy^-2" << endl;
    cout << "doseThres: " << doseThres << " Gy" << endl;
    cout << "arrestTime: " << arrestTime << " h" << endl;
    cout << "hypNecThres: " << hypNecThres <<" mmHg" << endl;
    cout << "dose: " << dose << " Gy" << endl;
    cout << "DO2: " << DO2 << " um^2/ms" << endl;
    cout << "Vmax: " << Vmax << " mmHg/ms" << endl;
    cout << "Km: " << Km << " mmHg" << endl;
    cout << "Dvegf: " << Dvegf << " um^2/ms" << endl;
    cout << "VmaxVegf: " << VmaxVegf << " mol/um^3ms" << endl;
    cout << "KmVegf: " << KmVegf << " mol/um^3" << endl;
    cout << "pO2NormVes: " << pO2NormVes << " mmHg" << endl;
    cout << "pO2TumVes: " << pO2TumVes << " mmHg" << endl;
    cout << "hypThres: " << hypThres << " mmHg" << endl;
    cout << "hypVegf: " << hypVegf << " mol/um^3" << endl;


    Treatment *treatment;
    treatment = new Treatment(dose, 80.0, 24.0, 0);

    Coupler *coupler;
    Tissue *model1;
    OxyTissue *model2;

    model1 = new Tissue(nrow, ncol, nlayer, inTum, inVes,
                        tumGrowth, tumTime, edgeOrder, cycDur, cycDistrib,
                        res, fibTime, ang, vascTumTime, vegfThres, alpha,
                        beta, doseThres, arrestTime, treatment, hypNecThres);

    double simTimeStep, oxySimTimeStep, sclFac, simTime;

    simTimeStep    = 6.0; //h;
    oxySimTimeStep = 1.0; //ms;
    sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
    
    double cellSize(20.0); //um
    double nu(oxySimTimeStep / (cellSize * cellSize));
    DO2 *= nu;
    Vmax *= oxySimTimeStep;
    Dvegf *= nu;
    VmaxVegf *= oxySimTimeStep;


    model2 = new OxyTissue(nrow, ncol, nlayer, inVes,
                           Dvegf, DO2, Vmax, Km,
                           pO2NormVes, pO2TumVes,
                           hypThres, VmaxVegf,
                           KmVegf, hypVegf);

    coupler = new Coupler(model1, model2);

    RootSimulator *sim;

    sim = new RootSimulator(coupler, simTimeStep,
                            oxySimTimeStep, sclFac);
    simTime = treatment->getDuration();
    sim->initSim();
    sim->simulate(simTimeStep, simTime);
    sim->stop();

    double finTumDens(model1->getOut()->at(0));
    double timeTo95(model1->getOut()->at(12));
    double timeTo99(model1->getOut()->at(13));
    double intTumDens(model1->getOut()->at(22));

    delete treatment;
    delete model1;
    delete model2;
    delete coupler;
    delete sim;

    cout << "finTumDens: " << finTumDens << endl;
    cout << "intTumDens: " << intTumDens << endl;
    cout << "timeTo95: " << timeTo95 << " h" << endl;
    cout << "timeTo99: " << timeTo99 << " h" << endl;

    vector<double> outputs(4);
    outputs[0] = finTumDens;
    outputs[1] = intTumDens;
    outputs[2] = timeTo95;
    outputs[3] = timeTo99;

    return outputs;
}


void morris(int K, int p, int N){
    int M(K + 1);
    vector<vector<double> > B(M, vector<double>(K, 0.0));

    for(int m(1); m < M; m++){
        for(int k(0); k <m ; k++){
            B[m][k] = 1.0;
        }
    }

    vector<int> vP;

    for(int k(0); k < K; k++){
        vP.push_back(k);
    }

    ifstream frefParInt("../InputFiles/refParInt.dat");
    double _pm1(1.0 / (p - 1.0));
    double delta(0.5 * _pm1 * p);
    double _delta(1.0 / delta);
    vector<double> h(K), x0(K);

    for(int k(0); k < K; k++){
        frefParInt >> x0[k];
        frefParInt >> h[k];
        h[k] -= x0[k];
    }
    frefParInt.close();

    int diffk, nEv(0), nEvTot((K + 1) * N);
    vector<double> xp(M);
    vector<vector<double> > Bp;
    vector<vector<double> > Bpp(M, vector<double>(K, 0.0));
    vector<vector<double> > y(M, vector<double>(4));

    vector<vector<double> > elEffTumDens(K, vector<double>(N, 0.0));
    vector<vector<double> > elEffIntTumDens(K, vector<double>(N, 0.0));
    vector<vector<double> > elEffTimeTo95(K, vector<double>(N, 0.0));
    vector<vector<double> > elEffTimeTo99(K, vector<double>(N, 0.0));

    vector<double> meanTumDens(K, 0.0);
    vector<double> meanIntTumDens(K, 0.0);
    vector<double> meanTimeTo95(K, 0.0);
    vector<double> meanTimeTo99(K, 0.0);

    ofstream fTumDens("../OutputFiles/tumDens.res");
    ofstream fIntTumDens("../OutputFiles/intTumDens.res");
    ofstream fTimeTo95("../OutputFiles/timeTo95.res");
    ofstream fTimeTo99("../OutputFiles/timeTo99.res");

    for(int n(0); n < N; n++){
        for(int m(0); m < M; m++){
            xp[m] = _pm1 * (rand() % ((p - 2) / 2 + 1));
        }

        bool perm;
        Bp = B;
        for(int k(0); k < K; k++){
            perm = rand() % 2;
            for(int m(0); m < M; m++){
                if(perm){
                    Bp[m][k] = 1.0 - Bp[m][k];
                }
                Bp[m][k] = xp[k] + delta * Bp[m][k];
            }
        }

        random_shuffle(vP.begin(), vP.end());

        for(int m(0); m < M; m++){
            for(int k(0); k < K; k++){
                Bpp[m][k] = x0[k] + Bp[m][vP[k]] * h[k];
            }
        }

        for(int m(0); m < M; m++){
            y[m] = model(Bpp[m]);
            nEv++;

            fTumDens    << y[m][0] << endl;
            fIntTumDens << y[m][1] << endl;
            fTimeTo95   << y[m][2] << endl;
            fTimeTo99   << y[m][3] << endl;

            cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
            cout << "---------------------------------------------" << endl;
        }

        for(int m(1); m < M; m++){
            diffk = find(vP.begin(), vP.end(), m - 1) - vP.begin();

            elEffTumDens[diffk][n]    = fabs(_delta * (y[m][0] - y[m - 1][0]));
            elEffIntTumDens[diffk][n] = fabs(_delta * (y[m][1] - y[m - 1][1]));
            elEffTimeTo95[diffk][n]   = fabs(_delta * (y[m][2] - y[m - 1][2]));
            elEffTimeTo99[diffk][n]   = fabs(_delta * (y[m][3] - y[m - 1][3]));

            meanTumDens[diffk]    += elEffTumDens[diffk][n];
            meanIntTumDens[diffk] += elEffIntTumDens[diffk][n];
            meanTimeTo95[diffk]   += elEffTimeTo95[diffk][n];
            meanTimeTo99[diffk]   += elEffTimeTo99[diffk][n];
        }
    }

    fTumDens.close();
    fIntTumDens.close();
    fTimeTo99.close();
    fTimeTo95.close();

    for(int k(0); k < K; k++){
        meanTumDens[k]    /= N;
        meanIntTumDens[k] /= N;
        meanTimeTo95[k]   /= N;
        meanTimeTo99[k]   /= N;
    }

    vector<double> stdDevTumDens(K, 0.0);
    vector<double> stdDevIntTumDens(K, 0.0);
    vector<double> stdDevTimeTo95(K, 0.0);
    vector<double> stdDevTimeTo99(K, 0.0);

    for(int k(0); k < K; k++){
        for(int n(0); n < N; n++){
            stdDevTumDens[k] += (elEffTumDens[k][n] - meanTumDens[k]) *
                    (elEffTumDens[k][n] - meanTumDens[k]);
            stdDevIntTumDens[k] += (elEffIntTumDens[k][n] - meanIntTumDens[k]) *
                    (elEffIntTumDens[k][n] - meanIntTumDens[k]);
            stdDevTimeTo95[k] += (elEffTimeTo95[k][n] - meanTimeTo95[k]) *
                    (elEffTimeTo95[k][n] - meanTimeTo95[k]);
            stdDevTimeTo99[k] += (elEffTimeTo99[k][n] - meanTimeTo99[k]) *
                    (elEffTimeTo99[k][n] - meanTimeTo99[k]);

        }

        stdDevTumDens[k]    = sqrt(stdDevTumDens[k] / (N - 1.0));
        stdDevIntTumDens[k] = sqrt(stdDevIntTumDens[k] / (N - 1.0));
        stdDevTimeTo95[k]   = sqrt(stdDevTimeTo95[k] / (N - 1.0));
        stdDevTimeTo99[k]   = sqrt(stdDevTimeTo99[k] / (N - 1.0));
    }

    ofstream fMorrisTumDens("../OutputFiles/morrisTumDens.res");
    ofstream fMorrisIntTumDens("../OutputFiles/morrisIntTumDens.res");
    ofstream fMorrisTimeTo95("../OutputFiles/morrisTimeTo95.res");
    ofstream fMorrisTimeTo99("../OutputFiles/morrisTimeTo99.res");

    for(int k(0); k < K; k++){
        fMorrisTumDens    << meanTumDens[k]    << " " << stdDevTumDens[k]    << endl;
        fMorrisIntTumDens << meanIntTumDens[k] << " " << stdDevIntTumDens[k] << endl;
        fMorrisTimeTo95   << meanTimeTo95[k]   << " " << stdDevTimeTo95[k]   << endl;
        fMorrisTimeTo99   << meanTimeTo99[k]   << " " << stdDevTimeTo99[k]   << endl;

        cout << meanTumDens[k]    << " " << stdDevTumDens[k]    << endl;
        cout << meanIntTumDens[k] << " " << stdDevIntTumDens[k] << endl;
        cout << meanTimeTo95[k]   << " " << stdDevTimeTo95[k]   << endl;
        cout << meanTimeTo99[k]   << " " << stdDevTimeTo99[k]   << endl;
    }
    fMorrisTumDens.close();
    fMorrisIntTumDens.close();
    fMorrisTimeTo95.close();
    fMorrisTimeTo99.close();
}


void morrisk(int kp, int K, int p, int N){
    ifstream frefParMean("../InputFiles/refParMean.dat");
    double _pm1(1.0 / (p - 1.0));
    double delta(0.5 * _pm1 * p);
    double _delta(1.0 / delta);
    double h, x0;
    vector<double> xmean(K);

    for(int k(0); k < K; k++){
        frefParMean >> xmean[k];
    }
    frefParMean.close();

    ifstream frefParInt("../InputFiles/refParInt.dat");
    for(int k(-1); k < kp; k++){
        frefParInt >> x0;
        frefParInt >> h;
    }
    frefParInt.close();

    h -= x0;

    int nEv(0), nEvTot(2 * N);
    vector<double> xp(2);
    vector<double> B(2);
    vector<vector<double> > Bp(2, vector<double>(K, 0.0));
    vector<vector<double> > y(2, vector<double>(4));

    vector<double> elEffTumDens(N);
    vector<double> elEffIntTumDens(N);
    vector<double> elEffTimeTo95(N);
    vector<double> elEffTimeTo99(N);

    double meanTumDens(0.0), meanIntTumDens(0.0);
    double meanTimeTo95(0.0), meanTimeTo99(0.0);

    ofstream fTumDens("../OutputFiles/tumDens.res");
    ofstream fIntTumDens("../OutputFiles/intTumDens.res");
    ofstream fTimeTo95("../OutputFiles/timeTo95.res");
    ofstream fTimeTo99("../OutputFiles/timeTo99.res");

    for(int n(0); n < N; n++){
        xp[0] = _pm1 * (rand() % ((p - 2) / 2 + 1));
        xp[1] = _pm1 * (rand() % ((p - 2) / 2 + 1));

        B[0] = xp[0];
        B[1] = xp[1] + delta;

        for(int m(0); m < 2; m++){
            for(int k(0); k < K; k++){
                Bp[m][k] = xmean[k];
            }
        }
        Bp[0][kp] = x0 + B[0] * h;
        Bp[1][kp] = x0 + B[1] * h;

        for(int m(0); m < 2; m++){
            y[m] = model(Bp[m]);
            nEv++;

            fTumDens    << y[m][0] << endl;
            fIntTumDens << y[m][1] << endl;
            fTimeTo95   << y[m][2] << endl;
            fTimeTo99   << y[m][3] << endl;

            cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
            cout << "---------------------------------------------" << endl;
        }

        elEffTumDens[n]    = fabs(_delta * (y[1][0] - y[0][0]));
        elEffIntTumDens[n] = fabs(_delta * (y[1][1] - y[0][1]));
        elEffTimeTo95[n]   = fabs(_delta * (y[1][2] - y[0][2]));
        elEffTimeTo99[n]   = fabs(_delta * (y[1][3] - y[0][3]));

        meanTumDens    += elEffTumDens[n];
        meanIntTumDens += elEffIntTumDens[n];
        meanTimeTo95   += elEffTimeTo95[n];
        meanTimeTo99   += elEffTimeTo99[n];
    }

    fTumDens.close();
    fIntTumDens.close();
    fTimeTo99.close();
    fTimeTo95.close();

    meanTumDens    /= N;
    meanIntTumDens /= N;
    meanTimeTo95   /= N;
    meanTimeTo99   /= N;

    double stdDevTumDens(0.0), stdDevIntTumDens(0.0);
    double stdDevTimeTo95(0.0), stdDevTimeTo99(0.0);

    for(int n(0); n < N; n++){
        stdDevTumDens += (elEffTumDens[n] - meanTumDens) *
                (elEffTumDens[n] - meanTumDens);
        stdDevIntTumDens += (elEffIntTumDens[n] - meanIntTumDens) *
                (elEffIntTumDens[n] - meanIntTumDens);
        stdDevTimeTo95 += (elEffTimeTo95[n] - meanTimeTo95) *
                (elEffTimeTo95[n] - meanTimeTo95);
        stdDevTimeTo99 += (elEffTimeTo99[n] - meanTimeTo99) *
                (elEffTimeTo99[n] - meanTimeTo99);

    }

    stdDevTumDens    = sqrt(stdDevTumDens / (N - 1.0));
    stdDevIntTumDens = sqrt(stdDevIntTumDens / (N - 1.0));
    stdDevTimeTo95   = sqrt(stdDevTimeTo95 / (N - 1.0));
    stdDevTimeTo99   = sqrt(stdDevTimeTo99 / (N - 1.0));

    ofstream fMorrisTumDens("../OutputFiles/morrisTumDens.res");
    ofstream fMorrisIntTumDens("../OutputFiles/morrisIntTumDens.res");
    ofstream fMorrisTimeTo95("../OutputFiles/morrisTimeTo95.res");
    ofstream fMorrisTimeTo99("../OutputFiles/morrisTimeTo99.res");

    fMorrisTumDens    << meanTumDens    << " " << stdDevTumDens    << endl;
    fMorrisIntTumDens << meanIntTumDens << " " << stdDevIntTumDens << endl;
    fMorrisTimeTo95   << meanTimeTo95   << " " << stdDevTimeTo95   << endl;
    fMorrisTimeTo99   << meanTimeTo99   << " " << stdDevTimeTo99   << endl;

    cout << meanTumDens    << " " << stdDevTumDens    << endl;
    cout << meanIntTumDens << " " << stdDevIntTumDens << endl;
    cout << meanTimeTo95   << " " << stdDevTimeTo95   << endl;
    cout << meanTimeTo99   << " " << stdDevTimeTo99   << endl;

    fMorrisTumDens.close();
    fMorrisIntTumDens.close();
    fMorrisTimeTo95.close();
    fMorrisTimeTo99.close();
}


void sobol(int K, int N){
    int nEv(0), nEvTot((K + 2) * N);
    vector<vector<double> > Xa(N, vector<double>(K));
    vector<vector<double> > Xb(N, vector<double>(K));
    vector<vector<double> > Xc(N, vector<double>(K));
    vector<double> Ya(N), Yb(N), Yc(N);

    for(int i(0); i < N; i++){
        //Matrices for a model considering the Legendre polynomial of degree d
        //Xa[i][0] = rand() % 5 + 1;
        //Xa[i][1] = -1.0 + 2.0 * (double)(rand()) / (double)(RAND_MAX);
        //Xb[i][0] = rand() % 5 + 1;
        //Xb[i][1] = -1.0 + 2.0 * (double)(rand()) / (double)(RAND_MAX);
        for(int j(0); j < K; j++){
            Xa[i][j] = (double)(rand()) / (double)(RAND_MAX);
            Xb[i][j] = (double)(rand()) / (double)(RAND_MAX);
        }
    }

    for(int i(0); i < N; i++){
        Ya[i] = toyModel(Xa[i]);
        nEv++;
        //cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
        //cout << "---------------------------------------------" << endl;
        Yb[i] = toyModel(Xb[i]);
        nEv++;
        //cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
        //cout << "---------------------------------------------" << endl;
    }

    double alpha, beta, sigma2, f0;
    vector<double> SI(K), TSI(K);
    ofstream fSens("../OutputFiles/sobol.res");

    fSens << K << endl;

    for(int k(0); k < K; k++){
        Xc = Xa;
        alpha = 0.0;
        beta = 0.0;
        sigma2 = 0.0;
        f0 = 0.0;
        for(int i(0); i < N; i++){
            Xc[i][k] = Xb[i][k];
            Yc[i] = toyModel(Xc[i]);
            nEv++;
            //cout << nEv << " out of " << nEvTot << " evaluations of the model";
            //cout << "---------------------------------------------" << endl;
            alpha += Yb[i] * Yc[i];
            beta += (Ya[i] - Yc[i]) * (Ya[i] - Yc[i]);
            f0 += Yb[i] + Yc[i];
            sigma2 += Yb[i] * Yb[i] + Yc[i] * Yc[i];
        }
        alpha /= N;
        f0 /= 2.0 * N;
        sigma2 /= 2.0 * N;
        sigma2 -= f0 * f0;
        SI[k] = (alpha - f0 * f0) / sigma2;
        TSI[k] = beta / (2.0 * N * sigma2);
        fSens << SI[k] << " " << TSI[k] << endl;
        cout << SI[k] << " " << TSI[k] << endl;
    }

    fSens.close();
}


void sobolFromFiles(int K){
    int N;
    vector<double> Ya, Yb, Yc;

    ifstream fYa("../OutputSobol/Ya.txt");
    double temp;
    fYa >> temp;
    while(!fYa.eof()){
        Ya.push_back(temp);
        fYa >> temp;
    }
    fYa.close();

    N = Ya.size();

    ifstream fYb("../OutputSobol/Yb.txt");
    fYb >> temp;
    while(!fYb.eof()){
        Yb.push_back(temp);
        fYb >> temp;
    }
    fYb.close();

    double alpha, beta, sigma2, f0;
    vector<double> SI(K), TSI(K);
    ofstream fSens("../OutputFiles/sobol.res");

    fSens << K << endl;

    for(int k(0); k < K; k++){
        ifstream fYc("../OutputSobol/Yc" + to_string(k) + ".txt");
        Yc.resize(N);
        alpha = 0.0;
        beta = 0.0;
        sigma2 = 0.0;
        f0 = 0.0;
        for(int i(0); i < N; i++){
            fYc >> Yc[i];
            alpha += Yb[i] * Yc[i];
            beta += (Ya[i] - Yc[i]) * (Ya[i] - Yc[i]);
            f0 += Yb[i] + Yc[i];
            sigma2 += Yb[i] * Yb[i] + Yc[i] * Yc[i];
        }
        alpha /= N;
        f0 /= 2.0 * N;
        sigma2 /= 2.0 * N;
        sigma2 -= f0 * f0;
        SI[k] = (alpha - f0 * f0) / sigma2;
        TSI[k] = beta / (2.0 * N * sigma2);
        fSens << SI[k] << " " << TSI[k] << endl;
        fYc.close();
    }

    fSens.close();
}


double toyModel(vector<double> x){
    //Legendre polynomial of degree d
    /*int d(x[0]);

    switch (d){
    case 1:
        return x[1];

    case 2:
        return 0.5 * (3.0 * x[1] * x[1] - 1.0);

    case 3:
        return 0.5 * (5.0 * pow(x[1], 3) -
                3.0 * x[1]);

    case 4:
        return 0.125 * (35.0 * pow(x[1], 4) -
                30.0 * x[1] * x[1] + 3.0);

    case 5:
        return 0.125 * (63.0 * pow(x[1], 5) -
                70.0 * pow(x[1], 3) +
                15.0 * x[1]);
    }*/
    //return x[0] * x[1];
    return x[0] + x[1]  + x[2] * x[2] + x[3] * x[4];
}
