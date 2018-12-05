#include "sensAn.hpp"

using namespace std;

void sobol(const int K, const int N, const int nOut,
           const double *x0, const double *h,
           double **SI, double **TSI,
           double ***SIConv, double ***TSIConv){
    int nEv(0), nEvTot((K + 2) * N);
    double **Xa, **Xb, **Xc;
    double **Ya, **Yb, **Yc;

    Xa = alloc2D(N, K);
    Xb = alloc2D(N, K);
    Xc = alloc2D(N, K);

    Ya = alloc2D(N, nOut);
    Yb = alloc2D(N, nOut);
    Yc = alloc2D(N, nOut);

    for(int i(0); i < N; i++){
        //Matrices for a model considering the Legendre polynomial of degree d
        //Xa[i][0] = rand() % 5 + 1;
        //Xa[i][1] = -1.0 + 2.0 * (double)(rand()) / (double)(RAND_MAX);
        //Xb[i][0] = rand() % 5 + 1;
        //Xb[i][1] = -1.0 + 2.0 * (double)(rand()) / (double)(RAND_MAX);
        for(int k(0); k < K; k++){
            Xa[i][k] = x0[k] + (double)(rand()) / (double)(RAND_MAX) * h[k];
            Xb[i][k] = x0[k] + (double)(rand()) / (double)(RAND_MAX) * h[k];
        }
    }

    for(int i(0); i < N; i++){
        toyModel(Xa[i], Ya[i]);
        //model(Xa[i], Ya[i]);
        nEv++;
        //cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
        //cout << "---------------------------------------------" << endl;
        toyModel(Xb[i], Yb[i]);
        //model(Xb[i], Yb[i]);
        nEv++;
        //cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
        //cout << "---------------------------------------------" << endl;
    }

    int iConv(0), nConv;
    double alpha[nOut], beta[nOut], sigma2[nOut], f0[nOut];
    double f02, f02Conv, sigma2Conv;

    for(int i(0); i < N; i++){
        for(int k(0); k < K; k++){
            Xc[i][k] = Xa[i][k];
        }
    }

    double _N(1.0 / N);

    nConv = 2;

    for(int j(0); j < nOut; j++){
        alpha[j]  = 0.0;
        beta[j]   = 0.0;
        f0[j]     = 0.0;
        sigma2[j] = 0.0;
    }

    for(int i(0); i < N; i++){
        Xc[i][0] = Xb[i][0];

        toyModel(Xc[i], Yc[i]);
        //model(Xc[i], Yc[i]);
        nEv++;
        //cout << nEv << " out of " << nEvTot << " evaluations of the model";
        //cout << "---------------------------------------------" << endl;

        for(int j(0); j < nOut; j++){
            alpha[j]  += Yb[i][j] * Yc[i][j];
            beta[j]   += (Ya[i][j] - Yc[i][j]) * (Ya[i][j] - Yc[i][j]);
            f0[j]     += Yb[i][j] + Yc[i][j];
            sigma2[j] += Yb[i][j] * Yb[i][j] + Yc[i][j] * Yc[i][j];
        }

        if(i == nConv){
            for(int j(0); j < nOut; j++){
                f02Conv    = 0.25 / (nConv * nConv) * f0[j] * f0[j] ;
                sigma2Conv = 0.5 / nConv * sigma2[j] - f02Conv;
                SIConv[iConv][j][0]  = (alpha[j] / nConv - f02Conv) / sigma2Conv;
                TSIConv[iConv][j][0] = 0.5 * beta[j] / (nConv * sigma2Conv);
            }
            iConv++;
            nConv *= 2;
        }
    }

    for(int j(0); j < nOut; j++){
        f02 = 0.25 * _N * _N * f0[j] * f0[j];
        sigma2[j] = 0.5 * _N * sigma2[j] - f02;
        SI[j][0] = (_N * alpha[j] - f02) / sigma2[j];
        TSI[j][0] = 0.5 * _N *  beta[j] / sigma2[j];
    }

    for(int k(1); k < K; k++){
        iConv = 0;
        nConv = 2;

        for(int j(0); j < nOut; j++){
            alpha[j]  = 0.0;
            beta[j]   = 0.0;
            f0[j]     = 0.0;
            sigma2[j] = 0.0;
        }

        for(int i(0); i < N; i++){
            Xc[i][k - 1] = Xa[i][k - 1];
            Xc[i][k] = Xb[i][k];

            toyModel(Xc[i], Yc[i]);
            //model(Xc[i], Yc[i]);
            nEv++;
            //cout << nEv << " out of " << nEvTot << " evaluations of the model";
            //cout << "---------------------------------------------" << endl;

            for(int j(0); j < nOut; j++){
                alpha[j]  += Yb[i][j] * Yc[i][j];
                beta[j]   += (Ya[i][j] - Yc[i][j]) * (Ya[i][j] - Yc[i][j]);
                f0[j]     += Yb[i][j] + Yc[i][j];
                sigma2[j] += Yb[i][j] * Yb[i][j] + Yc[i][j] * Yc[i][j];
            }

            if(i == nConv){
                for(int j(0); j < nOut; j++){
                    f02Conv    = 0.25 / (nConv * nConv) * f0[j] * f0[j] ;
                    sigma2Conv = 0.5 / nConv * sigma2[j] - f02Conv;
                    SIConv[iConv][j][k]  = (alpha[j] / nConv - f02Conv) / sigma2Conv;
                    TSIConv[iConv][j][k] = 0.5 * beta[j] / (nConv * sigma2Conv);
                }
                iConv++;
                nConv *= 2;
            }
        }

        for(int j(0); j < nOut; j++){
            f02       = 0.25 * _N * _N * f0[j] * f0[j] ;
            sigma2[j] = 0.5 * _N * sigma2[j] - f02;
            SI[j][k]  = (_N * alpha[j] - f02) / sigma2[j];
            TSI[j][k] = 0.5 * _N *  beta[j] / sigma2[j];
        }
    }

    free2D(Xa, N);
    free2D(Xb, N);
    free2D(Xc, N);

    free2D(Ya, N);
    free2D(Yb, N);
    free2D(Yc, N);
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

void sobolRT(const int N, const string nFRefParInt){
    const int K(34), NConv(log(N) / log(2.0)), nOut(6);
    double h[K], x0[K];
    ifstream fRefParInt(nFRefParInt.c_str());

    for(int k(0); k < K; k++){
        fRefParInt >> x0[k];
        fRefParInt >> h[k];
        h[k] -= x0[k];
    }

    fRefParInt.close();

    double **SI, **TSI;
    double ***SIConv, ***TSIConv;

    SI  = alloc2D(nOut, K);
    TSI = alloc2D(nOut, K);
    SIConv  = alloc3D(NConv, nOut, K);
    TSIConv = alloc3D(NConv, nOut, K);

    sobol(K, N, nOut, x0, h, SI, TSI, SIConv, TSIConv);

    ofstream fSobolEndTreatTumDens("../OutputFiles/sobolEndTreatTumDens.res");
    ofstream fSobol3MonTumDens("../OutputFiles/sobol3MonTumDens.res");
    ofstream fSobolRecTumDens("../OutputFiles/sobolRecTumDens.res");
    ofstream fSobolTumVol("../OutputFiles/sobolTumVol.res");
    ofstream fSobolIntTumDens("../OutputFiles/sobolIntTumDens.res");
    ofstream fSobolTimeTo95("../OutputFiles/sobolTimeTo95.res");
    ofstream fSobolTimeTo99("../OutputFiles/sobolTimeTo99.res");
    ofstream fSobolRecTime("../OutputFiles/sobolRecTime.res");

    ofstream fConvSIEndTreatTumDens("../OutputFiles/convSIEndTreatTumDens.res");
    ofstream fConvTSIEndTreatTumDens("../OutputFiles/convTSIEndTreatTumDens.res");
    ofstream fConvSI3MonTumDens("../OutputFiles/convSI3MonTumDens.res");
    ofstream fConvTSI3MonTumDens("../OutputFiles/convTS3MonTumDens.res");
    ofstream fConvSIRecTumDens("../OutputFiles/convSIRecTumDens.res");
    ofstream fConvTSIRecTumDens("../OutputFiles/convTSIRecTumDens.res");
    ofstream fConvSITumVol("../OutputFiles/convSITumVol.res");
    ofstream fConvTSITumVol("../OutputFiles/convTSITumVol.res");
    ofstream fConvSIIntTumDens("../OutputFiles/convSIIntTumDens.res");
    ofstream fConvTSIIntTumDens("../OutputFiles/convTSIIntTumDens.res");
    ofstream fConvSITimeTo95("../OutputFiles/convSITimeTo95.res");
    ofstream fConvTSITimeTo95("../OutputFiles/convTSITimeTo95.res");
    ofstream fConvSITimeTo99("../OutputFiles/convSITimeTo99.res");
    ofstream fConvTSITimeTo99("../OutputFiles/convTSITimeTo99.res");
    ofstream fConvSIRecTime("../OutputFiles/convSIRecTime.res");
    ofstream fConvTSIRecTime("../OutputFiles/convTSIRecTime.res");

    fSobolEndTreatTumDens << K << " " << 0 << endl;
    fSobol3MonTumDens     << K << " " << 0 << endl;
    fSobolRecTumDens      << K << " " << 0 << endl;
    fSobolTumVol          << K << " " << 0 << endl;
    fSobolIntTumDens      << K << " " << 0 << endl;
    fSobolTimeTo95        << K << " " << 0 << endl;
    fSobolTimeTo99        << K << " " << 0 << endl;
    fSobolRecTime         << K << " " << 0 << endl;

    for(int k(0); k < K; k++){
        fSobolEndTreatTumDens << SI[0][k] << " " << TSI[0][k] << endl;
        fSobol3MonTumDens     << SI[1][k] << " " << TSI[1][k] << endl;
        fSobolRecTumDens      << SI[2][k] << " " << TSI[2][k] << endl;
        fSobolTumVol          << SI[3][k] << " " << TSI[3][k] << endl;
        fSobolIntTumDens      << SI[4][k] << " " << TSI[4][k] << endl;
        fSobolTimeTo95        << SI[5][k] << " " << TSI[5][k] << endl;
        fSobolTimeTo99        << SI[6][k] << " " << TSI[6][k] << endl;
        fSobolRecTime         << SI[7][k] << " " << TSI[7][k] << endl;

        for(int l(0); l < NConv; l++){
            fConvSIEndTreatTumDens  << SIConv[l][0][k]  << " ";
            fConvTSIEndTreatTumDens << TSIConv[l][0][k] << " ";
            fConvSI3MonTumDens      << SIConv[l][1][k]  << " ";
            fConvTSI3MonTumDens     << TSIConv[l][1][k] << " ";
            fConvSIRecTumDens       << SIConv[l][2][k]  << " ";
            fConvTSIRecTumDens      << TSIConv[l][2][k] << " ";
            fConvSITumVol           << SIConv[l][3][k]  << " ";
            fConvTSITumVol          << TSIConv[l][3][k] << " ";
            fConvSIIntTumDens       << SIConv[l][4][k]  << " ";
            fConvTSIIntTumDens      << TSIConv[l][4][k] << " ";
            fConvSITimeTo95         << SIConv[l][5][k]  << " ";
            fConvTSITimeTo95        << TSIConv[l][5][k] << " ";
            fConvSITimeTo99         << SIConv[l][6][k]  << " ";
            fConvTSITimeTo99        << TSIConv[l][6][k] << " ";
            fConvSIRecTime          << SIConv[l][7][k]  << " ";
            fConvTSIRecTime         << TSIConv[l][7][k] << " ";
        }

        fConvSIEndTreatTumDens  << endl;
        fConvTSIEndTreatTumDens << endl;
        fConvSI3MonTumDens      << endl;
        fConvTSI3MonTumDens     << endl;
        fConvSIRecTumDens       << endl;
        fConvTSIRecTumDens      << endl;
        fConvSITumVol           << endl;
        fConvTSITumVol          << endl;
        fConvSIIntTumDens       << endl;
        fConvTSIIntTumDens      << endl;
        fConvSITimeTo95         << endl;
        fConvTSITimeTo95        << endl;
        fConvSITimeTo99         << endl;
        fConvTSITimeTo99        << endl;
        fConvSIRecTime          << endl;
        fConvTSIRecTime         << endl;
    }

    fSobolEndTreatTumDens.close();
    fSobol3MonTumDens.close();
    fSobolRecTumDens.close();
    fSobolTumVol.close();
    fSobolIntTumDens.close();
    fSobolTimeTo95.close();
    fSobolTimeTo99.close();
    fSobolRecTime.close();

    fConvSIEndTreatTumDens.close();
    fConvTSIEndTreatTumDens.close();
    fConvSI3MonTumDens.close();
    fConvTSI3MonTumDens.close();
    fConvSIRecTumDens.close();
    fConvTSIRecTumDens.close();
    fConvSITumVol.close();
    fConvTSITumVol.close();
    fConvSIIntTumDens.close();
    fConvTSIIntTumDens.close();
    fConvSITimeTo95.close();
    fConvTSITimeTo95.close();
    fConvSITimeTo99.close();
    fConvTSITimeTo99.close();
    fConvSIRecTime.close();
    fConvTSIRecTime.close();

    free2D(SI, nOut);
    free2D(TSI, nOut);
    free3D(SIConv, NConv, nOut);
    free3D(SIConv, NConv, nOut);
}


void sobolToy(const int N, const string nFRefParInt){
    const int K(5), NConv(log(N) / log(2.0)), nOut(1);
    double h[K], x0[K];
    ifstream fRefParInt(nFRefParInt.c_str());

    for(int k(0); k < K; k++){
        fRefParInt >> x0[k];
        fRefParInt >> h[k];
        h[k] -= x0[k];
    }

    fRefParInt.close();

    double **SI, **TSI;
    double ***SIConv, ***TSIConv;

    SI  = alloc2D(nOut, K);
    TSI = alloc2D(nOut, K);
    SIConv  = alloc3D(NConv, nOut, K);
    TSIConv = alloc3D(NConv, nOut, K);

    sobol(K, N, nOut, x0, h,  SI, TSI, SIConv, TSIConv);

    ofstream fSobolY("../OutputFiles/sobolY.res");
    ofstream fConvSIY("../OutputFiles/convSIY.res");
    ofstream fConvTSIY("../OutputFiles/convTSIY.res");

    fSobolY << K << " " << 0 << endl;

    for(int k(0); k < K; k++){
        fSobolY << SI[0][k] << " " << TSI[0][k] << endl;

        for(int l(0); l < NConv; l++){
            fConvSIY  << SIConv[l][0][k]  << " ";
            fConvTSIY << TSIConv[l][0][k] << " ";
        }

        fConvSIY  << endl;
        fConvTSIY << endl;
    }

    fSobolY.close();
    fConvSIY.close();
    fConvTSIY.close();

    free2D(SI, nOut);
    free2D(TSI, nOut);
    free3D(SIConv, NConv, nOut);
    free3D(TSIConv, NConv, nOut);
}
