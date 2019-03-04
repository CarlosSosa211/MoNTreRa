#include "sensAn.hpp"

using namespace std;

void morris(const int K, const int L, const int N, const int nOut, const int p,
            const double *x0, const double *h, double **mu, double **sigma,
            const string nFInTissueDim, const string nFInTum,
            const string nFInVes){
    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    if(!nFInTissueDim.empty() && !nFInTum.empty() && !nFInVes.empty()){
        readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                    cellSize, inTum, inVes);
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
            }
        }

        for(int m(0); m < M; m++){
            toyModel(Bp[m], y[m]);
            //model(Bp[m], y[m], nrow, ncol, nlayer, cellSize, inTum, inVes);
            nEv++;
            cout << nEv << " out of " << nEvTot <<
                    " evaluations of the model" << endl;
            cout << "---------------------------------------------" << endl;
        }

        for(int m(1); m < M; m++){
            diffk = find(vP.begin(), vP.end(), m - 1) - vP.begin();
            for(int j(0); j < nOut; j++){
                elEff[j][diffk][n] = fabs(_delta * (y[m][j] - y[m - 1][j]));
                mu[j][diffk] += elEff[j][diffk][n];
            }
        }
    }

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


void morrisRT(const int N, const int p, const string nFRefParInt,
              const string nFInTissueDim, const string nFInTum,
              const string nFInVes){
    const int K(35), nOut(15);
    //const int K(20), nOut(15);
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


void morrisToy(const int p, const int N, const string nFRefParInt){
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
