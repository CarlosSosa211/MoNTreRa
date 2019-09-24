#include "morris.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * This functions performs a Morris analysis.
 *
 * Inputs:
 *  - K: number of parameters,
 *  - L: number of repetitions of the analysis (number of times that the
 *  function is called; it is used to calculate the total number of evaluations,
 *  - N: number of repetitions,
 *  - nOut: number of outputs of the model,
 *  - p: number of levels of the Morris analysis,
 *  - x0: array containing the inferior values of the ranges of the parameters,
 *  - h: array containing the length of the ranges of the parameters,
 *  - nFInTissueDim: name of the file containing the dimensions of a tissue,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration of a tissue,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration of a tissue.
 *
 * Outputs:
 *  - mu: matrix containing the obtained mu* values; each row corresponds to an
 *  output and each column, to a parameter,
 *  - sigma: matrix containing the obtained sigma values; each row corresponds
 *  to an output and each column, to a parameter.
------------------------------------------------------------------------------*/

void morris(const int K, const int L, const int N, const int nOut, const int p,
            const double *x0, const double *h, double **mu, double **sigma,
            const string nFInTissueDim, const string nFInTum,
            const string nFInVes){
    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    if(!nFInTissueDim.empty() && !nFInTum.empty() && !nFInVes.empty()){
        readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol,
                    nlayer, cellSize, inTum, inVes);
    }
    else if(!nFInTissueDim.empty()){
        readInFiles(nFInTissueDim, nrow, ncol, nlayer, cellSize);
    }

    vector<int> vP;

    for(int k(0); k < K; k++){
        vP.push_back(k);
    }

    int diffk, nEv(0);
    const int M(K + 1);
    const int nEvTot(L * (K + 1) * N);
    const double _pm1(1.0 / (p - 1.0));
    const double delta(0.5 * _pm1 * p);
    const double _delta(1.0 / delta);
    double xp[M];
    double **B, **Bp;
    double **y;
    double ***elEff;
    ofstream fMorrisPar("../OutputFiles/morrisPar.res");
    ofstream fMorrisVarPar("../OutputFiles/morrisVarPar.res");
    ofstream fMorrisOut("../OutputFiles/morrisOut.res");

    B     = alloc2D(M, K);
    Bp    = alloc2D(M, K);
    y     = alloc2D(M, nOut);
    elEff = alloc3D(nOut, K, N);

    bool perm;

    for(int n(0); n < N; n++){
        for(int m(0); m < M; m++){
            xp[m] = _pm1 * (rand() % ((p - 2) / 2 + 1));
            for(int k(0); k < K; k++){
                B[m][k] = 0.0;
            }
        }

        for(int m(1); m < M; m++){
            for(int k(0); k < m ; k++){
                B[m][k] = 1.0;
            }
        }

        for(int k(0); k < K; k++){
            perm = rand() % 2;
            for(int m(0); m < M; m++){
                if(perm){
                    B[m][k] = 1.0 - B[m][k];
                }
                B[m][k] = xp[k] + delta * B[m][k];
            }
        }

        random_shuffle(vP.begin(), vP.end());

        for(int m(0); m < M; m++){
            for(int k(0); k < K; k++){
                Bp[m][k] = x0[k] + B[m][vP[k]] * h[k];
                fMorrisPar << Bp[m][k] << " ";
            }
            fMorrisPar << endl;
        }

        for(int m(0); m < M; m++){
            //toyModel(Bp[m], y[m]);
            //model(Bp[m], y[m], nrow, ncol, nlayer, cellSize, inTum, inVes);
            model(Bp[m], y[m], nrow, ncol, nlayer, cellSize);
            nEv++;
            cout << nEv << " out of " << nEvTot <<
                    " evaluations of the model" << endl;
            cout << "---------------------------------------------" << endl;
            for(int j(0); j < nOut; j++){
                fMorrisOut << y[m][j] << " ";
            }
            fMorrisOut << endl;
        }

        for(int m(1); m < M; m++){
            diffk = find(vP.begin(), vP.end(), m - 1) - vP.begin();
            fMorrisVarPar << diffk << endl;
            for(int j(0); j < nOut; j++){
                elEff[j][diffk][n] = fabs(_delta * (y[m][j] - y[m - 1][j]));
                mu[j][diffk] += elEff[j][diffk][n];
            }
        }
    }
    fMorrisVarPar.close();
    fMorrisPar.close();
    fMorrisOut.close();

    free2D(B, M);
    free2D(Bp, M);
    free2D(y, M);

    for(int k(0); k < K; k++){
        for(int j(0); j < nOut; j++){
            mu[j][k] /= N;
        }
    }

    for(int k(0); k < K; k++){
        for(int n(0); n < N; n++){
            for(int j(0); j < nOut; j++){
                sigma[j][k] += (elEff[j][k][n] - mu[j][k]) *
                        (elEff[j][k][n] - mu[j][k]);
            }
        }

        for(int j(0); j < nOut; j++){
            sigma[j][k] = sqrt(sigma[j][k] / (N - 1.0));
        }
    }

    free3D(elEff, nOut, K);
}


/*------------------------------------------------------------------------------
 * This functions performs a Morris analysis of the model of tumour
 * growth and response to radiotherapy and write in files the obtained results.
 *
 * Inputs:
 *  - N: number of repetitions.
------------------------------------------------------------------------------*/
void morrisFromFiles(const int N, const int p){
    const int K(42), M(K + 1), nOut(15);
    const double _pm1(1.0 / (p - 1.0));
    const double delta(0.5 * _pm1 * p);
    const double _delta(1.0 / delta);
    double ***X, ***Y;

    X = alloc3D(N, M, K);
    Y = alloc3D(N, M, nOut);

    ifstream fX("../OutputFiles/morrisPar.res");
    ifstream fY("../OutputFiles/morrisOut.res");

    for(int n(0); n < N; n++){
        for(int m(0); m < M; m++){
            for(int k(0); k < K; k++){
                fX >> X[n][m][k];
            }
            for(int j(0); j < nOut; j++){
                fY >> Y[n][m][j];
            }
        }
    }
    fX.close();
    fY.close();

    double ***elEff;

    elEff = alloc3D(nOut, K, N);

    double **mu, **sigma;

    mu    = alloc2D(nOut, K);
    sigma = alloc2D(nOut, K);
    for(int j(0); j < nOut; j++){
        for(int k(0); k < K; k++){
            mu[j][k]    = 0.0;
            sigma[j][k] = 0.0;
        }
    }
    bool diff(false);
    int diffk(0);

    int constk[] = {5, 7, 9, 34};
    int countConstk(0);
    for(int n(0); n < N; n++){
        countConstk = 0;
        for(int m(1); m < M; m++){
            diffk = 0;
            diff = X[n][m][diffk] != X[n][m - 1][diffk];
            while(diff == false && diffk < K){
                diffk++;
                diff = X[n][m][diffk] != X[n][m - 1][diffk];
            }
            if(!diff){
                diffk = constk[countConstk];
                countConstk++;

            }
            for(int j(0); j < nOut; j++){
                elEff[j][diffk][n] = fabs(_delta * (Y[n][m][j] -
                                                    Y[n][m - 1][j]));
                mu[j][diffk] += elEff[j][diffk][n];
            }
        }
    }
    free3D(X, N, M);
    free3D(Y, N, M);
    for(int k(0); k < K; k++){
        for(int j(0); j < nOut; j++){
            mu[j][k] /= N;
        }
    }
    for(int k(0); k < K; k++){
        for(int n(0); n < N; n++){
            for(int j(0); j < nOut; j++){
                sigma[j][k] += (elEff[j][k][n] - mu[j][k]) *
                        (elEff[j][k][n] - mu[j][k]);
            }
        }

        for(int j(0); j < nOut; j++){
            sigma[j][k] = sqrt(sigma[j][k] / (N - 1.0));
        }
    }

    free3D(elEff, nOut, K);

    ofstream fMorrisEndTreatTumDens("../OutputFiles/morrisEndTreatTumDens.res");
    ofstream fMorris3MonTumDens("../OutputFiles/morris3MonTumDens.res");
    ofstream fMorrisTumVol("../OutputFiles/morrisTumVol.res");
    ofstream fMorrisIntTumDens("../OutputFiles/morrisIntTumDens.res");
    ofstream fMorrisKilled50("../OutputFiles/morrisKilled50.res");
    ofstream fMorrisKilled80("../OutputFiles/morrisKilled80.res");
    ofstream fMorrisKilled90("../OutputFiles/morrisKilled90.res");
    ofstream fMorrisKilled95("../OutputFiles/morrisKilled95.res");
    ofstream fMorrisTimeTo95("../OutputFiles/morrisTimeTo95.res");
    ofstream fMorrisKilled99("../OutputFiles/morrisKilled99.res");
    ofstream fMorrisTimeTo99("../OutputFiles/morrisTimeTo99.res");
    ofstream fMorrisKilled999("../OutputFiles/morrisKilled999.res");
    ofstream fMorrisRec("../OutputFiles/morrisRec.res");
    ofstream fMorrisRecTumDens("../OutputFiles/morrisRecTumDens.res");
    ofstream fMorrisRecTime("../OutputFiles/morrisRecTime.res");

    for(int k(0); k < K; k++){
        fMorrisEndTreatTumDens << mu[0][k]  << " " << sigma[0][k]  << endl;
        fMorris3MonTumDens     << mu[1][k]  << " " << sigma[1][k]  << endl;
        fMorrisTumVol          << mu[2][k]  << " " << sigma[2][k]  << endl;
        fMorrisIntTumDens      << mu[3][k]  << " " << sigma[3][k]  << endl;
        fMorrisKilled50        << mu[4][k]  << " " << sigma[4][k]  << endl;
        fMorrisKilled80        << mu[5][k]  << " " << sigma[5][k]  << endl;
        fMorrisKilled90        << mu[6][k]  << " " << sigma[6][k]  << endl;
        fMorrisKilled95        << mu[7][k]  << " " << sigma[7][k]  << endl;
        fMorrisTimeTo95        << mu[8][k]  << " " << sigma[8][k]  << endl;
        fMorrisKilled99        << mu[9][k]  << " " << sigma[9][k]  << endl;
        fMorrisTimeTo99        << mu[10][k] << " " << sigma[10][k] << endl;
        fMorrisKilled999       << mu[11][k] << " " << sigma[11][k] << endl;
        fMorrisRec             << mu[12][k] << " " << sigma[12][k] << endl;
        fMorrisRecTumDens      << mu[13][k] << " " << sigma[13][k] << endl;
        fMorrisRecTime         << mu[14][k] << " " << sigma[14][k] << endl;
    }

    fMorrisEndTreatTumDens.close();
    fMorris3MonTumDens.close();
    fMorrisTumVol.close();
    fMorrisIntTumDens.close();
    fMorrisKilled50.close();
    fMorrisKilled80.close();
    fMorrisKilled90.close();
    fMorrisKilled95.close();
    fMorrisTimeTo95.close();
    fMorrisKilled99.close();
    fMorrisTimeTo99.close();
    fMorrisKilled999.close();
    fMorrisRec.close();
    fMorrisRecTumDens.close();
    fMorrisRecTime.close();

    free2D(mu, nOut);
    free2D(sigma, nOut);

}


/*------------------------------------------------------------------------------
 * This functions prepares and performs a Morris analysis of the model of tumour
 * growth and response to radiotherapy and write in files the obtained results.
 *
 * Inputs:
 *  - N: number of repetitions,
 *  - p: number of levels of the Morris analysis,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration of a tissue; by default, it is empty,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration of a tissue; by default, it is empty.
------------------------------------------------------------------------------*/

void morrisRT(const int N, const int p, const string nFRefParInt,
              const string nFInTissueDim, const string nFInTum,
              const string nFInVes){
    const int K(42), nOut(15);
    double h[K], x0[K];
    ifstream fRefParInt(nFRefParInt.c_str());

    for(int k(0); k < K; k++){
        fRefParInt >> x0[k];
        fRefParInt >> h[k];
        h[k] -= x0[k];
    }

    fRefParInt.close();

    double **mu, **sigma;

    mu    = alloc2D(nOut, K);
    sigma = alloc2D(nOut, K);

    for(int j(0); j < nOut; j++){
        for(int k(0); k < K; k++){
            mu[j][k]    = 0.0;
            sigma[j][k] = 0.0;
        }
    }

    morris(K, 1, N, nOut, p, x0, h, mu, sigma, nFInTissueDim, nFInTum, nFInVes);

    ofstream fMorrisEndTreatTumDens("../OutputFiles/morrisEndTreatTumDens.res");
    ofstream fMorris3MonTumDens("../OutputFiles/morris3MonTumDens.res");
    ofstream fMorrisTumVol("../OutputFiles/morrisTumVol.res");
    ofstream fMorrisIntTumDens("../OutputFiles/morrisIntTumDens.res");
    ofstream fMorrisKilled50("../OutputFiles/morrisKilled50.res");
    ofstream fMorrisKilled80("../OutputFiles/morrisKilled80.res");
    ofstream fMorrisKilled90("../OutputFiles/morrisKilled90.res");
    ofstream fMorrisKilled95("../OutputFiles/morrisKilled95.res");
    ofstream fMorrisTimeTo95("../OutputFiles/morrisTimeTo95.res");
    ofstream fMorrisKilled99("../OutputFiles/morrisKilled99.res");
    ofstream fMorrisTimeTo99("../OutputFiles/morrisTimeTo99.res");
    ofstream fMorrisKilled999("../OutputFiles/morrisKilled999.res");
    ofstream fMorrisRec("../OutputFiles/morrisRec.res");
    ofstream fMorrisRecTumDens("../OutputFiles/morrisRecTumDens.res");
    ofstream fMorrisRecTime("../OutputFiles/morrisRecTime.res");

    for(int k(0); k < K; k++){
        fMorrisEndTreatTumDens << mu[0][k]  << " " << sigma[0][k]  << endl;
        fMorris3MonTumDens     << mu[1][k]  << " " << sigma[1][k]  << endl;
        fMorrisTumVol          << mu[2][k]  << " " << sigma[2][k]  << endl;
        fMorrisIntTumDens      << mu[3][k]  << " " << sigma[3][k]  << endl;
        fMorrisKilled50        << mu[4][k]  << " " << sigma[4][k]  << endl;
        fMorrisKilled80        << mu[5][k]  << " " << sigma[5][k]  << endl;
        fMorrisKilled90        << mu[6][k]  << " " << sigma[6][k]  << endl;
        fMorrisKilled95        << mu[7][k]  << " " << sigma[7][k]  << endl;
        fMorrisTimeTo95        << mu[8][k]  << " " << sigma[8][k]  << endl;
        fMorrisKilled99        << mu[9][k]  << " " << sigma[9][k]  << endl;
        fMorrisTimeTo99        << mu[10][k] << " " << sigma[10][k] << endl;
        fMorrisKilled999       << mu[11][k] << " " << sigma[11][k] << endl;
        fMorrisRec             << mu[12][k] << " " << sigma[12][k] << endl;
        fMorrisRecTumDens      << mu[13][k] << " " << sigma[13][k] << endl;
        fMorrisRecTime         << mu[14][k] << " " << sigma[14][k] << endl;
    }

    fMorrisEndTreatTumDens.close();
    fMorris3MonTumDens.close();
    fMorrisTumVol.close();
    fMorrisIntTumDens.close();
    fMorrisKilled50.close();
    fMorrisKilled80.close();
    fMorrisKilled90.close();
    fMorrisKilled95.close();
    fMorrisTimeTo95.close();
    fMorrisKilled99.close();
    fMorrisTimeTo99.close();
    fMorrisKilled999.close();
    fMorrisRec.close();
    fMorrisRecTumDens.close();
    fMorrisRecTime.close();

    free2D(mu, nOut);
    free2D(sigma, nOut);
}


/*------------------------------------------------------------------------------
 * This functions prepares and performs a Morris analysis of the toy model and
 * write in files the obtained results.
 *
 * Inputs:
 *  - N: number of repetitions,
 *  - p: number of levels of the Morris analysis,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters.
------------------------------------------------------------------------------*/

void morrisToy(const int N, const int p, const string nFRefParInt){
    const int K(5), nOut(1);
    double h[K], x0[K];
    ifstream fRefParInt(nFRefParInt.c_str());

    for(int k(0); k < K; k++){
        fRefParInt >> x0[k];
        fRefParInt >> h[k];
        h[k] -= x0[k];
    }

    fRefParInt.close();

    double **mu, **sigma;

    mu    = alloc2D(nOut, K);
    sigma = alloc2D(nOut, K);

    for(int j(0); j < nOut; j++){
        for(int k(0); k < K; k++){
            mu[j][k]    = 0.0;
            sigma[j][k] = 0.0;
        }
    }

    morris(K, 1, N, nOut, p, x0, h, mu, sigma);

    ofstream fMorrisY("../OutputFiles/morrisY.res");

    for(int k(0); k < K; k++){
        fMorrisY << mu[0][k] << " " << sigma[0][k] << endl;
    }

    fMorrisY.close();

    free2D(mu, nOut);
    free2D(sigma, nOut);
}


/*------------------------------------------------------------------------------
 * This functions studies the impact of varying the range of a parameter when
 * performing a Morris analysis. Regular increments are considered.
 *
 * Inputs:
 *  - K: number of parameters,
 *  - kp: parameter of the model to be incremented,
 *  - L: number of regular increments,
 *  - N: number of repetitions,
 *  - nOut: number of outputs of the model,
 *  - p: number of levels of the Morris analysis,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration of a tissue,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration of a tissue.
 *
 * Outputs:
 *  - mu: 3D matrix containing the obtained mu* values; each row corresponds to
 *  an output, each column, to a parameter and each layer, to a range of
 *  parameter kp,
 *  - sigma: 3D matrix containing the obtained sigma values; each row
 *  corresponds to an output, each column, to a parameter and each layer, to a
 *  range of parameter kp.
------------------------------------------------------------------------------*/

void morrisVarRange(const int K, const int kp, const int L, const int N,
                    const int nOut, const int p, const string nFRefParInt,
                    double ***mu, double ***sigma, const string nFInTissueDim,
                    const string nFInTum, const string nFInVes){
    double h[K], x0[K];
    ifstream fRefParInt(nFRefParInt.c_str());

    for(int k(0); k < K; k++){
        fRefParInt >> x0[k];
        fRefParInt >> h[k];
        h[k] -= x0[k];
    }
    fRefParInt.close();

    int nEv(0);
    const int nEvTot((K + 1) * N * L);
    const double h10(0.1 * h[kp]);
    ofstream fVarRange("../OutputFiles/morrisVarRange.res");

    fVarRange << kp << endl;

    for(int l(0); l < L; l++){
        fVarRange << h[kp] << endl;
        morris(K, L, N, nOut, p, x0, h, mu[l], sigma[l], nFInTissueDim,
               nFInTum, nFInVes);
        h[kp] += h10;
    }

    fVarRange.close();
}


/*------------------------------------------------------------------------------
 * This functions prepares and performs a study of the impact of varying ranges
 * when performing a Morris analysis of the model of tumour growth and response
 * to radiotherapy and write in files the obtained results.
 *
 * Inputs:
 *  - kp: parameter of the model to be incremented,
 *  - L: number of regular increments,
 *  - N: number of repetitions,
 *  - p: number of levels of the Morris analysis,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration of a tissue,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration of a tissue,
 *  - nFInPO2: name of the file containing the initial pO2 values.
------------------------------------------------------------------------------*/

void morrisVarRangeRT(const int kp, const int L, const int N, const int p,
                      const string nFRefParInt, const string nFInTissueDim,
                      const string nFInTum, const string nFInVes){
    const int K(35), nOut(15);
    double ***mu, ***sigma;

    mu    = alloc3D(L, nOut, K);
    sigma = alloc3D(L, nOut, K);

    for(int l(0); l < L; l++){
        for(int j(0); j < nOut; j++){
            for(int k(0); k < K; k++){
                mu[l][j][k]    = 0.0;
                sigma[l][j][k] = 0.0;
            }
        }
    }

    morrisVarRange(K, kp, L, N, nOut, p, nFRefParInt, mu, sigma, nFInTissueDim,
                   nFInTum, nFInVes);

    for(int l(0); l < L; l++){
        ofstream fMorrisEndTreatTumDens("../OutputFiles/morrisEndTreatTumDens_"
                                        + to_string(l) + ".res");
        ofstream fMorris3MonTumDens("../OutputFiles/morris3MonTumDens_" +
                                    to_string(l) + ".res");
        ofstream fMorrisTumVol("../OutputFiles/morrisTumVol_" + to_string(l) +
                               ".res");
        ofstream fMorrisIntTumDens("../OutputFiles/morrisIntTumDens_" +
                                   to_string(l) + ".res");
        ofstream fMorrisKilled50("../OutputFiles/morrisKilled50_" +
                                 to_string(l) + ".res");
        ofstream fMorrisKilled80("../OutputFiles/morrisKilled80_" +
                                 to_string(l) + ".res");
        ofstream fMorrisKilled90("../OutputFiles/morrisKilled90_" +
                                 to_string(l) + ".res");
        ofstream fMorrisKilled95("../OutputFiles/morrisKilled95_" +
                                 to_string(l) + ".res");
        ofstream fMorrisTimeTo95("../OutputFiles/morrisTimeTo95_" +
                                 to_string(l) + ".res");
        ofstream fMorrisKilled99("../OutputFiles/morrisKilled99_" +
                                 to_string(l) + ".res");
        ofstream fMorrisTimeTo99("../OutputFiles/morrisTimeTo99_" +
                                 to_string(l) + ".res");
        ofstream fMorrisKilled999("../OutputFiles/morrisKilled999_" +
                                  to_string(l) + ".res");
        ofstream fMorrisRec("../OutputFiles/morrisRec_" + to_string(l) +
                            ".res");
        ofstream fMorrisRecTumDens("../OutputFiles/morrisRecTumDens_" +
                                   to_string(l) + ".res");
        ofstream fMorrisRecTime("../OutputFiles/morrisRecTime_" + to_string(l) +
                                ".res");

        for(int k(0); k < K; k++){
            fMorrisEndTreatTumDens << mu[l][0][k]  << " " << sigma[l][0][k]  <<
                                                      endl;
            fMorris3MonTumDens     << mu[l][1][k]  << " " << sigma[l][1][k]  <<
                                                      endl;
            fMorrisTumVol          << mu[l][2][k]  << " " << sigma[l][2][k]  <<
                                                      endl;
            fMorrisIntTumDens      << mu[l][3][k]  << " " << sigma[l][3][k]  <<
                                                      endl;
            fMorrisKilled50        << mu[l][4][k]  << " " << sigma[l][4][k]  <<
                                                      endl;
            fMorrisKilled80        << mu[l][5][k]  << " " << sigma[l][5][k]  <<
                                                      endl;
            fMorrisKilled90        << mu[l][6][k]  << " " << sigma[l][6][k]  <<
                                                      endl;
            fMorrisKilled95        << mu[l][7][k]  << " " << sigma[l][7][k]  <<
                                                      endl;
            fMorrisTimeTo95        << mu[l][8][k]  << " " << sigma[l][8][k]  <<
                                                      endl;
            fMorrisKilled99        << mu[l][9][k]  << " " << sigma[l][9][k]  <<
                                                      endl;
            fMorrisTimeTo99        << mu[l][10][k] << " " << sigma[l][10][k] <<
                                                      endl;
            fMorrisKilled999       << mu[l][11][k] << " " << sigma[l][11][k] <<
                                                      endl;
            fMorrisRec             << mu[l][12][k] << " " << sigma[l][12][k] <<
                                                      endl;
            fMorrisRecTumDens      << mu[l][13][k] << " " << sigma[l][13][k] <<
                                                      endl;
            fMorrisRecTime         << mu[l][14][k] << " " << sigma[l][14][k] <<
                                                      endl;
        }

        fMorrisEndTreatTumDens.close();
        fMorris3MonTumDens.close();
        fMorrisTumVol.close();
        fMorrisIntTumDens.close();
        fMorrisKilled50.close();
        fMorrisKilled80.close();
        fMorrisKilled90.close();
        fMorrisKilled95.close();
        fMorrisTimeTo95.close();
        fMorrisKilled99.close();
        fMorrisTimeTo99.close();
        fMorrisKilled999.close();
        fMorrisRec.close();
        fMorrisRecTumDens.close();
        fMorrisRecTime.close();
    }

    free3D(mu, L, nOut);
    free3D(sigma, L, nOut);
}


/*------------------------------------------------------------------------------
 * This functions prepares and performs a study of the impact of varying ranges
 * when performing a Morris analysis of the toy model and write in files the
 * obtained results.
 *
 * Inputs:
 *  - kp: parameter of the model to be incremented,
 *  - L: number of regular increments,
 *  - N: number of repetitions,
 *  - p: number of levels of the Morris analysis,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters.
------------------------------------------------------------------------------*/

void morrisVarRangeToy(const int kp, const int L, const int N, const int p,
                       const string nFRefParInt){
    const int K(5), nOut(1);
    double ***mu, ***sigma;

    mu    = alloc3D(L, nOut, K);
    sigma = alloc3D(L, nOut, K);

    for(int l(0); l < L; l++){
        for(int j(0); j < nOut; j++){
            for(int k(0); k < K; k++){
                mu[l][j][k]    = 0.0;
                sigma[l][j][k] = 0.0;
            }
        }
    }

    morrisVarRange(K, kp, L, N, nOut, p, nFRefParInt, mu, sigma);

    for(int l(0); l < L; l++){
        ofstream fMorrisY("../OutputFiles/morrisY_" + to_string(l) + ".res");
        for(int k(0); k < K; k++){
            fMorrisY << mu[l][0][k] << " " << sigma[l][0][k] << endl;
        }
        fMorrisY.close();
    }

    free3D(mu, L, nOut);
    free3D(sigma, L, nOut);
}
