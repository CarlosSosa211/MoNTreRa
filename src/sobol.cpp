#include "sobol.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * This function performs a Sobol analysis.
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
        readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol,
                    nlayer, cellSize, inTum, inVes);
    }

    int nEv(0), nEvTot((K + 2) * N);
    double **Xa, **Xb, **Xc;
    double **Ya, **Yb, **Yc;
    ofstream fXa("../OutputFiles/Xa.res");
    ofstream fXb("../OutputFiles/Xb.res");
    ofstream fYa("../OutputFiles/Ya.res");
    ofstream fYb("../OutputFiles/Yb.res");

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
            fXa << Xa[i][k] << " ";
            fXb << Xb[i][k] << " ";
        }
        fXa << endl;
        fXb << endl;
    }
    fXa.close();
    fXb.close();

    /*ifstream fXa("../OutputFiles/Xa.res");
    ifstream fXb("../OutputFiles/Xb.res");
    ifstream fYa("../OutputFiles/Ya.res");
    ifstream fYb("../OutputFiles/Yb.res");

    Xa = alloc2D(N, K);
    Xb = alloc2D(N, K);
    Xc = alloc2D(N, K);

    Ya = alloc2D(N, nOut);
    Yb = alloc2D(N, nOut);
    Yc = alloc2D(N, nOut);

    for(int i(0); i < N; i++){
        for(int k(0); k < K; k++){
            fXa >> Xa[i][k];
            fXb >> Xb[i][k];
        }
    }
    fXa.close();
    fXb.close();

    for(int i(0); i < N; i++){
        for(int j(0); j < nOut; j++){
            fYa >> Ya[i][j];
            fYb >> Yb[i][j];
        }
    }

    fYa.close();
    fYb.close();*/

    for(int i(0); i < N; i++){
        //toyModel(Xa[i], Ya[i]);
        //model(Xa[i], Ya[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        //reducedModel(Xa[i], Ya[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        reducedFracModel(Xa[i], Ya[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        nEv++;
        cout << nEv << " out of " << nEvTot << " evaluations of the model" <<
                endl;
        cout << "---------------------------------------------" << endl;
        for(int j(0); j < nOut; j++){
            fYa << Ya[i][j] << " ";
        }
        fYa << endl;

        //toyModel(Xb[i], Yb[i]);
        //model(Xb[i], Yb[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        //reducedModel(Xb[i], Yb[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        reducedFracModel(Xa[i], Ya[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        nEv++;
        cout << nEv << " out of " << nEvTot << " evaluations of the model" <<
                endl;
        cout << "---------------------------------------------" << endl;
        for(int j(0); j < nOut; j++){
            fYb << Yb[i][j] << " ";
        }
        fYb << endl;
    }
    fYa.close();
    fYb.close();

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

    ofstream fXc("../OutputFiles/Xc0.res");
    ofstream fYc("../OutputFiles/Yc0.res");

    for(int i(0); i < N; i++){
        Xc[i][0] = Xb[i][0];
        for(int k(0); k < K; k++){
            fXc << Xc[i][k] << " ";
        }
        fXc << endl;
        //toyModel(Xc[i], Yc[i]);
        //model(Xc[i], Yc[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        //reducedModel(Xc[i], Yc[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        reducedFracModel(Xa[i], Ya[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
        nEv++;
        cout << nEv << " out of " << nEvTot << " evaluations of the model";
        cout << "---------------------------------------------" << endl;
        for(int j(0); j < nOut; j++){
            fYc << Yc[i][j] << " ";
        }
        fYc << endl;

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

    fXc.close();
    fYc.close();
    for(int k(1); k < K; k++){
        fXc.open("../OutputFiles/Xc" + to_string(k) + ".res");
        fYc.open("../OutputFiles/Yc" + to_string(k) + ".res");
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
            for(int k(0); k < K; k++){
                fXc << Xc[i][k] << " ";
            }
            fXc << endl;

            //toyModel(Xc[i], Yc[i]);
            //model(Xc[i], Yc[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
            /*reducedModel(Xc[i], Yc[i], nrow, ncol, nlayer, cellSize, inTum,
                         inVes);*/
            reducedFracModel(Xa[i], Ya[i], nrow, ncol, nlayer, cellSize, inTum, inVes);
            nEv++;
            cout << nEv << " out of " << nEvTot << " evaluations of the model";
            cout << "---------------------------------------------" << endl;
            for(int j(0); j < nOut; j++){
                fYc << Yc[i][j] << " ";
            }
            fYc << endl;

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
        fXc.close();
        fYc.close();
    }


    free2D(Xa, N);
    free2D(Xb, N);
    free2D(Xc, N);

    free2D(Ya, N);
    free2D(Yb, N);
    free2D(Yc, N);
}


/*------------------------------------------------------------------------------
 * This function performs a Sobol analysis using outputs already calculated and
 * read from files.
 *
 * Inputs:
 *  - K: number of parameters.
------------------------------------------------------------------------------*/

void sobolFromFiles(){
    int K, N, nOut;
    double **Ya, **Yb, **Yc;

    ifstream fSobolPar("../OutputFiles/sobolPar.res");
    fSobolPar >> K >> N >> nOut;
    fSobolPar.close();

    Ya = alloc2D(N, nOut);
    Yb = alloc2D(N, nOut);
    Yc = alloc2D(N, nOut);

    ifstream fYa("../OutputFiles/Ya.res");
    for(int i(0); i < N; i++){
        for(int j(0); j < nOut; j++){
            fYa >> Ya[i][j];
        }
    }
    fYa.close();

    ifstream fYb("../OutputFiles/Yb.res");
    for(int i(0); i < N; i++){
        for(int j(0); j < nOut; j++){
            fYb >> Yb[i][j];
        }
    }
    fYb.close();

    double alpha, beta, sigma2, f0;
    double **SI, **TSI;
    SI  = alloc2D(nOut, K);
    TSI = alloc2D(nOut, K);

    ofstream fSobol8wTumDens("../OutputFiles/sobol8wTumDens.res");
    ofstream fSobol12wTumDens("../OutputFiles/sobol12wTumDens.res");
    ofstream fSobol8wTumVol("../OutputFiles/sobol8wTumVol.res");
    ofstream fSobol12wTumVol("../OutputFiles/sobol12wTumVol.res");
    ofstream fSobol8wIntTumDens("../OutputFiles/sobol8wIntTumDens.res");
    ofstream fSobol12wIntTumDens("../OutputFiles/sobol12wIntTumDens.res");
    ofstream fSobol8wIntTumVol("../OutputFiles/sobol8wIntTumVol.res");
    ofstream fSobol12wIntTumVol("../OutputFiles/sobol12wIntTumVol.res");
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
    ofstream fSobolRecTumVol("../OutputFiles/sobolRecTumVol.res");
    ofstream fSobolRecTime("../OutputFiles/sobolRecTime.res");

    for(int k(0); k < K; k++){
        ifstream fYc("../OutputFiles/Yc" + to_string(k) + ".res");
        for(int i(0); i < N; i++){
            for(int j(0); j < nOut; j++){
                fYc >> Yc[i][j];
            }
        }
        fYc.close();

        for(int j(0); j < nOut; j++){
            alpha  = 0.0;
            beta   = 0.0;
            sigma2 = 0.0;
            f0     = 0.0;

            for(int i(0); i < N; i++){
                alpha  += Yb[i][j] * Yc[i][j];
                beta   += (Ya[i][j] - Yc[i][j]) * (Ya[i][j] - Yc[i][j]);
                f0     += Yb[i][j] + Yc[i][j];
                sigma2 += Yb[i][j] * Yb[i][j] + Yc[i][j] * Yc[i][j];
            }

            alpha  /= N;
            f0     /= 2.0 * N;
            sigma2 /= 2.0 * N;
            sigma2 -= f0 * f0;
            SI[j][k]  = (alpha - f0 * f0) / sigma2;
            TSI[j][k] = beta / (2.0 * N * sigma2);
        }

        fSobol8wTumDens     << SI[0][k]  << " " << TSI[0][k]  << endl;
        fSobol12wTumDens    << SI[1][k]  << " " << TSI[1][k]  << endl;
        fSobol8wTumVol      << SI[2][k]  << " " << TSI[2][k]  << endl;
        fSobol12wTumVol     << SI[3][k]  << " " << TSI[3][k]  << endl;
        fSobol8wIntTumDens  << SI[4][k]  << " " << TSI[4][k]  << endl;
        fSobol12wIntTumDens << SI[5][k]  << " " << TSI[5][k]  << endl;
        fSobol8wIntTumVol   << SI[6][k]  << " " << TSI[6][k]  << endl;
        fSobol12wIntTumVol  << SI[7][k]  << " " << TSI[7][k]  << endl;
        fSobolKilled50      << SI[8][k]  << " " << TSI[8][k]  << endl;
        fSobolKilled80      << SI[9][k]  << " " << TSI[9][k]  << endl;
        fSobolKilled90      << SI[10][k] << " " << TSI[10][k] << endl;
        fSobolKilled95      << SI[11][k] << " " << TSI[11][k] << endl;
        fSobolTimeTo95      << SI[12][k] << " " << TSI[12][k] << endl;
        fSobolKilled99      << SI[13][k] << " " << TSI[13][k] << endl;
        fSobolTimeTo99      << SI[14][k] << " " << TSI[14][k] << endl;
        fSobolKilled999     << SI[15][k] << " " << TSI[15][k] << endl;
        fSobolRec           << SI[16][k] << " " << TSI[16][k] << endl;
        fSobolRecTumDens    << SI[17][k] << " " << TSI[17][k] << endl;
        fSobolRecTumVol     << SI[18][k] << " " << TSI[18][k] << endl;
        fSobolRecTime       << SI[19][k] << " " << TSI[19][k] << endl;
        fYc.close();
    }

    fSobol8wTumDens.close();
    fSobol12wTumDens.close();
    fSobol8wTumVol.close();
    fSobol12wTumVol.close();
    fSobol8wIntTumDens.close();
    fSobol12wIntTumDens.close();
    fSobol8wIntTumVol.close();
    fSobol12wIntTumVol.close();
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
    fSobolRecTumVol.close();
    fSobolRecTime.close();
}

/*------------------------------------------------------------------------------
 * This function prepares and performs a Sobol analysis of the model of tumour
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
    const int K(2), NConv(log(N) / log(2.0)), nOut(20);
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

    ofstream fSobol8wTumDens("../OutputFiles/sobol8wTumDens.res");
    ofstream fSobol12wTumDens("../OutputFiles/sobol12wTumDens.res");
    ofstream fSobol8wTumVol("../OutputFiles/sobol8wTumVol.res");
    ofstream fSobol12wTumVol("../OutputFiles/sobol12wTumVol.res");
    ofstream fSobol8wIntTumDens("../OutputFiles/sobol8wIntTumDens.res");
    ofstream fSobol12wIntTumDens("../OutputFiles/sobol12wIntTumDens.res");
    ofstream fSobol8wIntTumVol("../OutputFiles/sobol8wIntTumVol.res");
    ofstream fSobol12wIntTumVol("../OutputFiles/sobol12wIntTumVol.res");
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
    ofstream fSobolRecTumVol("../OutputFiles/sobolRecTumVol.res");
    ofstream fSobolRecTime("../OutputFiles/sobolRecTime.res");

    ofstream fConvSI8wTumDens("../OutputFiles/convSI8wTumDens.res");
    ofstream fConvTSI8wTumDens("../OutputFiles/convTSI8wTumDens.res");
    ofstream fConvSI12wTumDens("../OutputFiles/convSI12wTumDens.res");
    ofstream fConvTSI12wTumDens("../OutputFiles/convTSI12wTumDens.res");
    ofstream fConvSI8wTumVol("../OutputFiles/convSI8wTumVol.res");
    ofstream fConvTSI8wTumVol("../OutputFiles/convTSI8wTumVol.res");
    ofstream fConvSI12wTumVol("../OutputFiles/convSI12wTumVol.res");
    ofstream fConvTSI12wTumVol("../OutputFiles/convTSI12wTumVol.res");
    ofstream fConvSI8wIntTumDens("../OutputFiles/convSI8wIntTumDens.res");
    ofstream fConvTSI8wIntTumDens("../OutputFiles/convTSI8wIntTumDens.res");
    ofstream fConvSI12wIntTumDens("../OutputFiles/convSI12wIntTumDens.res");
    ofstream fConvTSI12wIntTumDens("../OutputFiles/convTSI12wIntTumDens.res");
    ofstream fConvSI8wIntTumVol("../OutputFiles/convSI8wIntTumVol.res");
    ofstream fConvTSI8wIntTumVol("../OutputFiles/convTSI8wIntTumVol.res");
    ofstream fConvSI12wIntTumVol("../OutputFiles/convSI12wIntTumVol.res");
    ofstream fConvTSI12wIntTumVol("../OutputFiles/convTSI12wIntTumVol.res");
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
    ofstream fConvSIRecTumVol("../OutputFiles/convSIRecTumVol.res");
    ofstream fConvTSIRecTumVol("../OutputFiles/convTSIRecTumVol.res");
    ofstream fConvSIRecTime("../OutputFiles/convSIRecTime.res");
    ofstream fConvTSIRecTime("../OutputFiles/convTSIRecTime.res");


    fSobol8wTumDens     << K << " " << 0 << endl;
    fSobol12wTumDens    << K << " " << 0 << endl;
    fSobol8wTumVol      << K << " " << 0 << endl;
    fSobol12wTumVol     << K << " " << 0 << endl;
    fSobol8wIntTumDens  << K << " " << 0 << endl;
    fSobol12wIntTumDens << K << " " << 0 << endl;
    fSobol8wIntTumVol   << K << " " << 0 << endl;
    fSobol12wIntTumVol  << K << " " << 0 << endl;
    fSobolKilled50      << K << " " << 0 << endl;
    fSobolKilled80      << K << " " << 0 << endl;
    fSobolKilled90      << K << " " << 0 << endl;
    fSobolKilled95      << K << " " << 0 << endl;
    fSobolTimeTo95      << K << " " << 0 << endl;
    fSobolKilled99      << K << " " << 0 << endl;
    fSobolTimeTo99      << K << " " << 0 << endl;
    fSobolKilled999     << K << " " << 0 << endl;
    fSobolRec           << K << " " << 0 << endl;
    fSobolRecTumDens    << K << " " << 0 << endl;
    fSobolRecTumVol     << K << " " << 0 << endl;
    fSobolRecTime       << K << " " << 0 << endl;

    for(int k(0); k < K; k++){
        fSobol8wTumDens     << SI[0][k]  << " " << TSI[0][k]  << endl;
        fSobol12wTumDens    << SI[1][k]  << " " << TSI[1][k]  << endl;
        fSobol8wTumVol      << SI[2][k]  << " " << TSI[2][k]  << endl;
        fSobol12wTumVol     << SI[3][k]  << " " << TSI[3][k]  << endl;
        fSobol8wIntTumDens  << SI[4][k]  << " " << TSI[4][k]  << endl;
        fSobol12wIntTumDens << SI[5][k]  << " " << TSI[5][k]  << endl;
        fSobol8wIntTumVol   << SI[6][k]  << " " << TSI[6][k]  << endl;
        fSobol12wIntTumVol  << SI[7][k]  << " " << TSI[7][k]  << endl;
        fSobolKilled50      << SI[8][k]  << " " << TSI[8][k]  << endl;
        fSobolKilled80      << SI[9][k]  << " " << TSI[9][k]  << endl;
        fSobolKilled90      << SI[10][k] << " " << TSI[10][k] << endl;
        fSobolKilled95      << SI[11][k] << " " << TSI[11][k] << endl;
        fSobolTimeTo95      << SI[12][k] << " " << TSI[12][k] << endl;
        fSobolKilled99      << SI[13][k] << " " << TSI[13][k] << endl;
        fSobolTimeTo99      << SI[14][k] << " " << TSI[14][k] << endl;
        fSobolKilled999     << SI[15][k] << " " << TSI[15][k] << endl;
        fSobolRec           << SI[16][k] << " " << TSI[16][k] << endl;
        fSobolRecTumDens    << SI[17][k] << " " << TSI[17][k] << endl;
        fSobolRecTumVol     << SI[18][k] << " " << TSI[18][k] << endl;
        fSobolRecTime       << SI[19][k] << " " << TSI[19][k] << endl;

        for(int l(0); l < NConv; l++){
            fConvSI8wTumDens      << SIConv[l][0][k]  << " ";
            fConvTSI8wTumDens     << TSIConv[l][0][k] << " ";
            fConvSI12wTumDens     << SIConv[l][1][k]  << " ";
            fConvTSI12wTumDens    << TSIConv[l][1][k] << " ";
            fConvSI8wTumVol       << SIConv[l][2][k]  << " ";
            fConvTSI8wTumVol      << TSIConv[l][2][k] << " ";
            fConvSI12wTumVol      << SIConv[l][3][k]  << " ";
            fConvTSI12wTumVol     << TSIConv[l][3][k] << " ";
            fConvSI8wIntTumDens   << SIConv[l][4][k]  << " ";
            fConvTSI8wIntTumDens  << TSIConv[l][4][k] << " ";
            fConvSI12wIntTumDens  << SIConv[l][5][k]  << " ";
            fConvTSI12wIntTumDens << TSIConv[l][5][k] << " ";
            fConvSI8wIntTumVol    << SIConv[l][6][k]  << " ";
            fConvTSI8wIntTumVol   << TSIConv[l][6][k] << " ";
            fConvSI12wIntTumVol   << SIConv[l][7][k]  << " ";
            fConvTSI12wIntTumVol  << TSIConv[l][7][k] << " ";
            fConvSIKilled50       << SIConv[l][8][k]  << " ";
            fConvTSIKilled50      << TSIConv[l][8][k] << " ";
            fConvSIKilled80       << SIConv[l][9][k]  << " ";
            fConvTSIKilled80      << TSIConv[l][9][k] << " ";
            fConvSIKilled90       << SIConv[l][10][k]  << " ";
            fConvTSIKilled90      << TSIConv[l][10][k] << " ";
            fConvSIKilled95       << SIConv[l][11][k]  << " ";
            fConvTSIKilled95      << TSIConv[l][11][k] << " ";
            fConvSITimeTo95       << SIConv[l][12][k]  << " ";
            fConvTSITimeTo95      << TSIConv[l][12][k] << " ";
            fConvSIKilled99       << SIConv[l][13][k]  << " ";
            fConvTSIKilled99      << TSIConv[l][13][k] << " ";
            fConvSITimeTo99       << SIConv[l][14][k]  << " ";
            fConvTSITimeTo99      << TSIConv[l][14][k] << " ";
            fConvSIKilled999      << SIConv[l][15][k]  << " ";
            fConvTSIKilled999     << TSIConv[l][15][k] << " ";
            fConvSIRec            << SIConv[l][16][k]  << " ";
            fConvTSIRec           << TSIConv[l][16][k] << " ";
            fConvSIRecTumDens     << SIConv[l][17][k]  << " ";
            fConvTSIRecTumDens    << TSIConv[l][17][k] << " ";
            fConvSIRecTumVol      << SIConv[l][18][k]  << " ";
            fConvTSIRecTumVol     << TSIConv[l][18][k] << " ";
            fConvSIRecTime        << SIConv[l][19][k]  << " ";
            fConvTSIRecTime       << TSIConv[l][19][k] << " ";
        }

        fConvSI8wTumDens      << endl;
        fConvTSI8wTumDens     << endl;
        fConvSI12wTumDens     << endl;
        fConvTSI12wTumDens    << endl;
        fConvSI8wTumVol       << endl;
        fConvTSI8wTumVol      << endl;
        fConvSI12wTumVol      << endl;
        fConvTSI12wTumVol     << endl;
        fConvSI8wIntTumDens   << endl;
        fConvTSI8wIntTumDens  << endl;
        fConvSI12wIntTumDens  << endl;
        fConvTSI12wIntTumDens << endl;
        fConvSI8wIntTumVol    << endl;
        fConvTSI8wIntTumVol   << endl;
        fConvSI12wIntTumVol   << endl;
        fConvTSI12wIntTumVol  << endl;
        fConvSIKilled50       << endl;
        fConvTSIKilled50      << endl;
        fConvSIKilled80       << endl;
        fConvTSIKilled80      << endl;
        fConvSIKilled90       << endl;
        fConvTSIKilled90      << endl;
        fConvSIKilled95       << endl;
        fConvTSIKilled95      << endl;
        fConvSITimeTo95       << endl;
        fConvTSITimeTo95      << endl;
        fConvSIKilled99       << endl;
        fConvTSIKilled99      << endl;
        fConvSITimeTo99       << endl;
        fConvTSITimeTo99      << endl;
        fConvSIKilled999      << endl;
        fConvTSIKilled999     << endl;
        fConvSIRec            << endl;
        fConvTSIRec           << endl;
        fConvSIRecTumDens     << endl;
        fConvTSIRecTumDens    << endl;
        fConvSIRecTumVol      << endl;
        fConvTSIRecTumVol     << endl;
        fConvSIRecTime        << endl;
        fConvTSIRecTime       << endl;
    }

    fSobol8wTumDens.close();
    fSobol12wTumDens.close();
    fSobol8wTumVol.close();
    fSobol12wTumVol.close();
    fSobol8wIntTumDens.close();
    fSobol12wIntTumDens.close();
    fSobol8wIntTumVol.close();
    fSobol12wIntTumVol.close();
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
    fSobolRecTumVol.close();
    fSobolRecTime.close();

    fConvSI8wTumDens.close();
    fConvTSI8wTumDens.close();
    fConvSI12wTumDens.close();
    fConvTSI12wTumDens.close();
    fConvSI8wTumVol.close();
    fConvTSI8wTumVol.close();
    fConvSI12wTumVol.close();
    fConvTSI12wTumVol.close();
    fConvSI8wIntTumDens.close();
    fConvTSI8wIntTumDens.close();
    fConvSI12wIntTumDens.close();
    fConvTSI12wIntTumDens.close();
    fConvSI8wIntTumVol.close();
    fConvTSI8wIntTumVol.close();
    fConvSI12wIntTumVol.close();
    fConvTSI12wIntTumVol.close();
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
    fConvSIRecTumVol.close();
    fConvTSIRecTumVol.close();
    fConvSIRecTime.close();
    fConvTSIRecTime.close();

    free2D(SI, nOut);
    free2D(TSI, nOut);
    free3D(SIConv, NConv, nOut);
    free3D(SIConv, NConv, nOut);
}


/*------------------------------------------------------------------------------
 * This function prepares and performs a Sobol analysis of the toy model and
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
