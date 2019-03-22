#include "sensAn.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * This functions performs a Sobol analysis.
 *
 * Inputs:
 *  - K: number of parameters,
 *  - N: number of repetitions,
 *  - nOut: number of outputs of the model,
 *  - x0: array containing the inferior values of the ranges of the parameters,
 *  - h: array containing the length of the ranges of the parameters,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration.
 *
 * Outputs:
 *  - SI: matrix containing the obtained SI values; each row corresponds to an
 *  output and each column, to a parameter,
 *  - TSI: matrix containing the obtained TSI values; each row corresponds
 *  to an output and each column, to a parameter,
 *  - SIConv: 3D matrix containing the obtained SI values; each row corresponds
 *  to an output, each column, to a parameter and each layer, to a number of
 *  repetitions of the form 2^i < N,
 *  - TSIConv: 3D matrix containing the obtained TSI values; each row
 *  corresponds to an output, each column, to a parameter and each layer, to a
 *  number of repetitions of the form 2^i < N.
------------------------------------------------------------------------------*/

void sobol(const int K, const int N, const int nOut, const double *x0,
           const double *h, double **SI, double **TSI, double ***SIConv,
           double ***TSIConv, const string nFInTissueDim, const string nFInTum,
           const string nFInVes){
    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    if(!nFInTissueDim.empty() && !nFInTum.empty() && !nFInVes.empty()){
        readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                    cellSize, inTum, inVes);
    }

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
        for(int k(0); k < K; k++){
            Xa[i][k] = x0[k] + (double)(rand()) / (double)(RAND_MAX) * h[k];
            Xb[i][k] = x0[k] + (double)(rand()) / (double)(RAND_MAX) * h[k];
        }
    }

    for(int i(0); i < N; i++){
        //toyModel(Xa[i], Ya[i]);
        model(Xa[i], Ya[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        nEv++;
        cout << nEv << " out of " << nEvTot << " evaluations of the model" <<
                endl;
        cout << "---------------------------------------------" << endl;
        //toyModel(Xb[i], Yb[i]);
        model(Xb[i], Yb[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        nEv++;
        cout << nEv << " out of " << nEvTot << " evaluations of the model" <<
                endl;
        cout << "---------------------------------------------" << endl;
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

        //toyModel(Xc[i], Yc[i]);
        model(Xc[i], Yc[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        nEv++;
        cout << nEv << " out of " << nEvTot << " evaluations of the model";
        cout << "---------------------------------------------" << endl;

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
                SIConv[iConv][j][0]  = (alpha[j] / nConv - f02Conv) /
                        sigma2Conv;
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

            //toyModel(Xc[i], Yc[i]);
            model(Xc[i], Yc[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
            nEv++;
            cout << nEv << " out of " << nEvTot << " evaluations of the model";
            cout << "---------------------------------------------" << endl;

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
                    SIConv[iConv][j][k]  = (alpha[j] / nConv - f02Conv) /
                            sigma2Conv;
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

/*------------------------------------------------------------------------------
 * This functions performs a Sobol analysis using outputs already calculated and
 * read from files.
 *
 * Inputs:
 *  - K: number of parameters.
------------------------------------------------------------------------------*/

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

/*------------------------------------------------------------------------------
 * This functions prepares and perform a Sobol analysis of the model of tumour
 * growth and response to radiotherapy and writes the obtained results in files.
 *
 * Inputs:
 *  - N: number of repetitions,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration.
------------------------------------------------------------------------------*/

void sobolRT(const int N, const string nFRefParInt, const string nFInTissueDim,
             const string nFInTum, const string nFInVes){
    const int K(34), NConv(log(N) / log(2.0)), nOut(15);
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

    sobol(K, N, nOut, x0, h, SI, TSI, SIConv, TSIConv, nFInTissueDim, nFInTum,
          nFInVes);

    ofstream fSobolEndTreatTumDens("../OutputFiles/sobolEndTreatTumDens.res");
    ofstream fSobol3MonTumDens("../OutputFiles/sobol3MonTumDens.res");
    ofstream fSobolTumVol("../OutputFiles/sobolTumVol.res");
    ofstream fSobolIntTumDens("../OutputFiles/sobolIntTumDens.res");
    ofstream fSobolKilled50("../OutputFiles/sobolKilled50.res");
    ofstream fSobolKilled80("../OutputFiles/sobolKilled80.res");
    ofstream fSobolKilled90("../OutputFiles/sobolKilled90.res");
    ofstream fSobolKilled95("../OutputFiles/sobolKilled95.res");
    ofstream fSobolTimeTo95("../OutputFiles/sobolTimeTo95.res");
    ofstream fSobolKilled99("../OutputFiles/sobolKilled99.res");
    ofstream fSobolTimeTo99("../OutputFiles/sobolTimeTo99.res");
    ofstream fSobolKilled999("../OutputFiles/sobolKilled999.res");
    ofstream fSobolRec("../OutputFiles/sobolRec.res");
    ofstream fSobolRecTumDens("../OutputFiles/sobolRecTumDens.res");
    ofstream fSobolRecTime("../OutputFiles/sobolRecTime.res");

    ofstream fConvSIEndTreatTumDens("../OutputFiles/convSIEndTreatTumDens.res");
    ofstream fConvTSIEndTreatTumDens("../OutputFiles/convTSIEndTreatTumDens.res");
    ofstream fConvSI3MonTumDens("../OutputFiles/convSI3MonTumDens.res");
    ofstream fConvTSI3MonTumDens("../OutputFiles/convTS3MonTumDens.res");
    ofstream fConvSITumVol("../OutputFiles/convSITumVol.res");
    ofstream fConvTSITumVol("../OutputFiles/convTSITumVol.res");
    ofstream fConvSIIntTumDens("../OutputFiles/convSIIntTumDens.res");
    ofstream fConvTSIIntTumDens("../OutputFiles/convTSIIntTumDens.res");
    ofstream fConvSIKilled50("../OutputFiles/convSIKilled50.res");
    ofstream fConvTSIKilled50("../OutputFiles/convTSIKilled50.res");
    ofstream fConvSIKilled80("../OutputFiles/convSIKilled80.res");
    ofstream fConvTSIKilled80("../OutputFiles/convTSIKilled80.res");
    ofstream fConvSIKilled90("../OutputFiles/convSIKilled90.res");
    ofstream fConvTSIKilled90("../OutputFiles/convTSIKilled90.res");
    ofstream fConvSIKilled95("../OutputFiles/convSIKilled95.res");
    ofstream fConvTSIKilled95("../OutputFiles/convTSIKilled95.res");
    ofstream fConvSITimeTo95("../OutputFiles/convSITimeTo95.res");
    ofstream fConvTSITimeTo95("../OutputFiles/convTSITimeTo95.res");
    ofstream fConvSIKilled99("../OutputFiles/convSIKilled99.res");
    ofstream fConvTSIKilled99("../OutputFiles/convTSIKilled99.res");
    ofstream fConvSITimeTo99("../OutputFiles/convSITimeTo99.res");
    ofstream fConvTSITimeTo99("../OutputFiles/convTSITimeTo99.res");
    ofstream fConvSIKilled999("../OutputFiles/convSIKilled999.res");
    ofstream fConvTSIKilled999("../OutputFiles/convTSIKilled999.res");
    ofstream fConvSIRec("../OutputFiles/convSIRec.res");
    ofstream fConvTSIRec("../OutputFiles/convTSIRec.res");
    ofstream fConvSIRecTumDens("../OutputFiles/convSIRecTumDens.res");
    ofstream fConvTSIRecTumDens("../OutputFiles/convTSIRecTumDens.res");
    ofstream fConvSIRecTime("../OutputFiles/convSIRecTime.res");
    ofstream fConvTSIRecTime("../OutputFiles/convTSIRecTime.res");

    fSobolEndTreatTumDens << K << " " << 0 << endl;
    fSobol3MonTumDens     << K << " " << 0 << endl;
    fSobolTumVol          << K << " " << 0 << endl;
    fSobolIntTumDens      << K << " " << 0 << endl;
    fSobolKilled50        << K << " " << 0 << endl;
    fSobolKilled80        << K << " " << 0 << endl;
    fSobolKilled90        << K << " " << 0 << endl;
    fSobolKilled95        << K << " " << 0 << endl;
    fSobolTimeTo95        << K << " " << 0 << endl;
    fSobolKilled99        << K << " " << 0 << endl;
    fSobolTimeTo99        << K << " " << 0 << endl;
    fSobolKilled999       << K << " " << 0 << endl;
    fSobolRec             << K << " " << 0 << endl;
    fSobolRecTumDens      << K << " " << 0 << endl;
    fSobolRecTime         << K << " " << 0 << endl;

    for(int k(0); k < K; k++){
        fSobolEndTreatTumDens << SI[0][k]  << " " << TSI[0][k]  << endl;
        fSobol3MonTumDens     << SI[1][k]  << " " << TSI[1][k]  << endl;
        fSobolTumVol          << SI[2][k]  << " " << TSI[2][k]  << endl;
        fSobolIntTumDens      << SI[3][k]  << " " << TSI[3][k]  << endl;
        fSobolKilled50        << SI[4][k]  << " " << TSI[4][k]  << endl;
        fSobolKilled80        << SI[5][k]  << " " << TSI[5][k]  << endl;
        fSobolKilled90        << SI[6][k]  << " " << TSI[6][k]  << endl;
        fSobolKilled95        << SI[7][k]  << " " << TSI[7][k]  << endl;
        fSobolTimeTo95        << SI[8][k]  << " " << TSI[8][k]  << endl;
        fSobolKilled99        << SI[9][k]  << " " << TSI[9][k]  << endl;
        fSobolTimeTo99        << SI[10][k] << " " << TSI[10][k] << endl;
        fSobolKilled999       << SI[11][k] << " " << TSI[11][k] << endl;
        fSobolRec             << SI[12][k] << " " << TSI[12][k] << endl;
        fSobolRecTumDens      << SI[13][k] << " " << TSI[13][k] << endl;
        fSobolRecTime         << SI[14][k] << " " << TSI[14][k] << endl;

        for(int l(0); l < NConv; l++){
            fConvSIEndTreatTumDens  << SIConv[l][0][k]  << " ";
            fConvTSIEndTreatTumDens << TSIConv[l][0][k] << " ";
            fConvSI3MonTumDens      << SIConv[l][1][k]  << " ";
            fConvTSI3MonTumDens     << TSIConv[l][1][k] << " ";
            fConvSITumVol           << SIConv[l][2][k]  << " ";
            fConvTSITumVol          << TSIConv[l][2][k] << " ";
            fConvSIIntTumDens       << SIConv[l][3][k]  << " ";
            fConvTSIIntTumDens      << TSIConv[l][3][k] << " ";
            fConvSIKilled50         << SIConv[l][4][k]  << " ";
            fConvTSIKilled50        << TSIConv[l][4][k] << " ";
            fConvSIKilled80         << SIConv[l][5][k]  << " ";
            fConvTSIKilled80        << TSIConv[l][5][k] << " ";
            fConvSIKilled90         << SIConv[l][6][k]  << " ";
            fConvTSIKilled90        << TSIConv[l][6][k] << " ";
            fConvSIKilled95         << SIConv[l][7][k]  << " ";
            fConvTSIKilled95        << TSIConv[l][7][k] << " ";
            fConvSITimeTo95         << SIConv[l][8][k]  << " ";
            fConvTSITimeTo95        << TSIConv[l][8][k] << " ";
            fConvSIKilled99         << SIConv[l][9][k]  << " ";
            fConvTSIKilled99        << TSIConv[l][9][k] << " ";
            fConvSITimeTo99         << SIConv[l][10][k]  << " ";
            fConvTSITimeTo99        << TSIConv[l][10][k] << " ";
            fConvSIKilled999        << SIConv[l][11][k]  << " ";
            fConvTSIKilled999       << TSIConv[l][11][k] << " ";
            fConvSIRec              << SIConv[l][12][k]  << " ";
            fConvTSIRec             << TSIConv[l][12][k] << " ";
            fConvSIRecTumDens       << SIConv[l][13][k]  << " ";
            fConvTSIRecTumDens      << TSIConv[l][13][k] << " ";
            fConvSIRecTime          << SIConv[l][14][k]  << " ";
            fConvTSIRecTime         << TSIConv[l][14][k] << " ";
        }

        fConvSIEndTreatTumDens  << endl;
        fConvTSIEndTreatTumDens << endl;
        fConvSI3MonTumDens      << endl;
        fConvTSI3MonTumDens     << endl;
        fConvSITumVol           << endl;
        fConvTSITumVol          << endl;
        fConvSIIntTumDens       << endl;
        fConvTSIIntTumDens      << endl;
        fConvSIKilled50         << endl;
        fConvTSIKilled50        << endl;
        fConvSIKilled80         << endl;
        fConvTSIKilled80        << endl;
        fConvSIKilled90         << endl;
        fConvTSIKilled90        << endl;
        fConvSIKilled95         << endl;
        fConvTSIKilled95        << endl;
        fConvSITimeTo95         << endl;
        fConvTSITimeTo95        << endl;
        fConvSIKilled99         << endl;
        fConvTSIKilled99        << endl;
        fConvSITimeTo99         << endl;
        fConvTSITimeTo99        << endl;
        fConvSIKilled999        << endl;
        fConvTSIKilled999       << endl;
        fConvSIRec              << endl;
        fConvTSIRec             << endl;
        fConvSIRecTumDens       << endl;
        fConvTSIRecTumDens      << endl;
        fConvSIRecTime          << endl;
        fConvTSIRecTime         << endl;
    }

    fSobolEndTreatTumDens.close();
    fSobol3MonTumDens.close();
    fSobolTumVol.close();
    fSobolIntTumDens.close();
    fSobolKilled50.close();
    fSobolKilled80.close();
    fSobolKilled90.close();
    fSobolKilled95.close();
    fSobolTimeTo95.close();
    fSobolKilled99.close();
    fSobolTimeTo99.close();
    fSobolKilled999.close();
    fSobolRec.close();
    fSobolRecTumDens.close();
    fSobolRecTime.close();

    fConvSIEndTreatTumDens.close();
    fConvTSIEndTreatTumDens.close();
    fConvSI3MonTumDens.close();
    fConvTSI3MonTumDens.close();
    fConvSITumVol.close();
    fConvTSITumVol.close();
    fConvSIIntTumDens.close();
    fConvTSIIntTumDens.close();
    fConvSIKilled50.close();
    fConvTSIKilled50.close();
    fConvSIKilled80.close();
    fConvTSIKilled80.close();
    fConvSIKilled90.close();
    fConvTSIKilled90.close();
    fConvSIKilled95.close();
    fConvTSIKilled95.close();
    fConvSITimeTo95.close();
    fConvTSITimeTo95.close();
    fConvSIKilled99.close();
    fConvTSIKilled99.close();
    fConvSITimeTo99.close();
    fConvTSITimeTo99.close();
    fConvSIKilled999.close();
    fConvTSIKilled999.close();
    fConvSIRec.close();
    fConvTSIRec.close();
    fConvSIRecTumDens.close();
    fConvTSIRecTumDens.close();
    fConvSIRecTime.close();
    fConvTSIRecTime.close();

    free2D(SI, nOut);
    free2D(TSI, nOut);
    free3D(SIConv, NConv, nOut);
    free3D(SIConv, NConv, nOut);
}


/*------------------------------------------------------------------------------
 * This functions prepares and perform a Sobol analysis of the toy model and
 * writes the obtained results in files.
 *
 * Inputs:
 *  - N: number of repetitions,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters.
------------------------------------------------------------------------------*/

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
