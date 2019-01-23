#include "sensAn.hpp"

using namespace std;

void var1ParRange(const int kp, const int L, const string nRefParInt,
                  const string nFInTissueDim, const string nFInTum,
                  const string nFInVes){
    const int K(39);
    double h[K], x0[K], x[K];
    ifstream fRefParInt(nRefParInt.c_str());

    for(int k(0); k < K; k++){
        fRefParInt >> x0[k];
        fRefParInt >> h[k];
        h[k] -= x0[k];
        x[k] = x0[k] + 0.5 * h[k];
    }
    x[kp] = x0[kp];
    fRefParInt.close();

    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                cellSize, inTum, inVes);

    const int nOut(8);
    const double delta(h[kp] / (L - 1));
    double y[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream fEndTreatTumDens, f3MonTumDens, fFinTumVol, fIntTumDens;
    ofstream fKilled50, fKilled80, fKilled90, fKilled95, fKilled99, fKilled999;
    ofstream fTimeTo95, fTimeTo99;
    ofstream fRec, fRecTumDens, fRecTime;

    for(int i(0); i < L; i++){
        nFTumDens     = "../OutputFiles/tumDens_" + to_string(i) + ".res";
        nFTumVol      = "../OutputFiles/tumVol_" + to_string(i) + ".res";
        nFVascDens    = "../OutputFiles/vascDens_" + to_string(i) + ".res";
        nFKilledCells = "../OutputFiles/killedCells_" + to_string(i) + ".res";
        nFCycle       = "../OutputFiles/cycle_" + to_string(i) + ".res";
        nFHypDens     = "../OutputFiles/hypDens_" + to_string(i) + ".res";
        nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(i) + ".res";
        nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(i) + ".res";

        fEndTreatTumDens.open("../OutputFiles/endTreatTumDens_" + to_string(i) + ".res");
        f3MonTumDens.open("../OutputFiles/3MonTumDens_" + to_string(i) + ".res");
        fFinTumVol.open("../OutputFiles/finTumVol_" + to_string(i) + ".res");
        fIntTumDens.open("../OutputFiles/intTumDens_" + to_string(i) + ".res");
        fKilled50.open("../OutputFiles/killed50_" + to_string(i) + ".res");
        fKilled80.open("../OutputFiles/killed80_" + to_string(i) + ".res");
        fKilled90.open("../OutputFiles/killed90_" + to_string(i) + ".res");
        fKilled95.open("../OutputFiles/killed95_" + to_string(i) + ".res");
        fTimeTo95.open("../OutputFiles/timeTo95_" + to_string(i) + ".res");
        fKilled99.open("../OutputFiles/killed99_" + to_string(i) + ".res");
        fTimeTo99.open("../OutputFiles/timeTo99_" + to_string(i) + ".res");
        fKilled999.open("../OutputFiles/killed999_" + to_string(i) + ".res");
        fRec.open("../OutputFiles/rec_" + to_string(i) + ".res");
        fRecTumDens.open("../OutputFiles/recTumDens_" + to_string(i) + ".res");
        fRecTime.open("../OutputFiles/recTime_" + to_string(i) + ".res");

        model(x, y, nrow, ncol, nlayer, cellSize, inTum, inVes, nFTumDens,
              nFTumVol, nFVascDens, nFKilledCells, nFCycle, nFHypDens,
              nFPO2Stat, nFVegfStat);

        fEndTreatTumDens << y[0];
        f3MonTumDens     << y[1];
        fFinTumVol       << y[2];
        fIntTumDens      << y[3];
        fKilled50        << y[4];
        fKilled80        << y[5];
        fKilled90        << y[6];
        fKilled95        << y[7];
        fTimeTo95        << y[8];
        fKilled99        << y[9];
        fTimeTo99        << y[10];
        fKilled999       << y[11];
        fRec             << y[12];
        fRecTumDens      << y[13];
        fRecTime         << y[14];

        fEndTreatTumDens.close();
        f3MonTumDens.close();
        fFinTumVol.close();
        fIntTumDens.close();
        fKilled50.close();
        fKilled80.close();
        fKilled90.close();
        fKilled95.close();
        fTimeTo95.close();
        fKilled99.close();
        fTimeTo99.close();
        fKilled999.close();
        fRec.close();
        fRecTumDens.close();
        fRecTime.close();

        cout << i + 1 << " out of " << L << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
        x[kp] += delta;
    }
}


void varErr(const string nFVarPar, const string nFMostRelPar, const string nFLeastPar,
            const string nFInTissueDim, const string nFInTum, const string nFInVes,
            const int L, const int P){
    const int K(38), nOut(15);
    int nLeastRelPar, nMostRelPar, nVarPar;

    ifstream fMostRelPar(nFMostRelPar.c_str());

    fMostRelPar >> nMostRelPar;

    int Xi[nMostRelPar];
    double h, x[K], X[nMostRelPar][L];

    for(int k(0); k < nMostRelPar; k++){
        fMostRelPar >> Xi[k] >> X[k][0] >> X[k][5];
        h = X[k][5] - X[k][0];
        X[k][1] = X[k][0] + 0.2 * h;
        X[k][2] = X[k][0] + 0.4 * h;
        X[k][3] = X[k][0] + 0.6 * h;
        X[k][4] = X[k][0] + 0.8 * h;
    }
    fMostRelPar.close();

    ifstream fLeastRelPar(nFLeastPar.c_str());

    fLeastRelPar >> nLeastRelPar;

    int i;
    for(int k(0); k < nLeastRelPar; k++){
        fLeastRelPar >> i;
        fLeastRelPar >> x[i];
    }
    fLeastRelPar.close();

    ifstream fVarPar(nFVarPar.c_str());

    fVarPar >> nVarPar;

    int XVari[nVarPar];
    double XVar[nVarPar][2];

    for(int k(0); k < nVarPar; k++){
        fVarPar >> XVari[k] >> XVar[k][0] >> XVar[k][1];
    }
    fVarPar.close();

    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                cellSize, inTum, inVes);

    int nEv(0), nEvTot(2 * pow(L, nMostRelPar) * P);
    double y0[P][nOut], y0mean[nOut], y0std[nOut];
    double y1[P][nOut], y1mean[nOut], y1std[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream fEndTreatTumDens("../OutputFiles/endTreatTumDens.res");
    ofstream f3MonTumDens("../OutputFiles/3MonTumDens.res");
    ofstream fFinTumVol("../OutputFiles/finTumVol.res");
    ofstream fIntTumDens("../OutputFiles/intTumDens.res");
    ofstream fKilled50("../OutputFiles/killed50.res");
    ofstream fKilled80("../OutputFiles/killed80.res");
    ofstream fKilled90("../OutputFiles/killed90.res");
    ofstream fKilled95("../OutputFiles/killed95.res");
    ofstream fTimeTo95("../OutputFiles/timeTo95.res");
    ofstream fKilled99("../OutputFiles/killed99.res");
    ofstream fTimeTo99("../OutputFiles/timeTo99.res");
    ofstream fKilled999("../OutputFiles/killed999.res");
    ofstream fRec("../OutputFiles/rec.res");
    ofstream fRecTumDens("../OutputFiles/recTumDens.res");
    ofstream fRecTime("../OutputFiles/recTime.res");
    ofstream fCombPar("../OutputFiles/combPar.res");

    int count(0);

    for(int i1(0); i1 < L; i1++){
        x[Xi[0]] = X[0][i1];
        for(int i2(0); i2 < L; i2++){
            x[Xi[1]] = X[1][i2];
            for(int i3(0); i3 < L; i3++){
                x[Xi[2]] = X[2][i3];
                for(int k(0); k < nVarPar; k++){
                    x[XVari[k]] = XVar[k][0];
                }

                for(int j(0); j < nOut; j++){
                    y0mean[j] = 0.0;
                    y0std[j]  = 0.0;
                    y1mean[j] = 0.0;
                    y1std[j]  = 0.0;
                }

                for(int p(0); p < P; p++){
                    nFTumDens     = "../OutputFiles/tumDens_" + to_string(count) + "_" +
                            "0_" + to_string(p) + ".res";
                    nFTumVol      = "../OutputFiles/tumVol_" + to_string(count) + "_" +
                            "0_" + to_string(p) + ".res";
                    nFVascDens    = "../OutputFiles/vascDens_" + to_string(count) + "_" +
                            "0_" + to_string(p) + ".res";
                    nFKilledCells = "../OutputFiles/killedCells_" + to_string(count) + "_" +
                            "0_" + to_string(p) + ".res";
                    nFCycle       = "../OutputFiles/cycle_" + to_string(count) + "_" +
                            "0_" + to_string(p) + ".res";
                    nFHypDens     = "../OutputFiles/hypDens_" + to_string(count) + "_" +
                            "0_" + to_string(p) + ".res";
                    nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(count) + "_" +
                            "0_" + to_string(p) + ".res";
                    nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(count) + "_" +
                            "0_" + to_string(p) + ".res";

                    model(x, y0[p],  nrow, ncol, nlayer, cellSize, inTum, inVes, nFTumDens,
                          nFTumVol, nFVascDens, nFKilledCells, nFCycle, nFHypDens,
                          nFPO2Stat, nFVegfStat);
                    nEv++;

                    cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
                    cout << "---------------------------------------------" << endl;

                    for(int j(0); j < nOut; j++){
                        y0mean[j] += y0[p][j];
                    }
                }

                for(int j(0); j < nOut; j++){
                    y0mean[j] /= P;
                }

                for(int k(0); k < nVarPar; k++){
                    x[XVari[k]] = XVar[k][1];
                }

                for(int p(0); p < P; p++){
                    nFTumDens     = "../OutputFiles/tumDens_" + to_string(count) + "_" +
                            "1_" + to_string(p) + ".res";
                    nFTumVol      = "../OutputFiles/tumVol_" + to_string(count) + "_" +
                            "1_" + to_string(p) + ".res";
                    nFVascDens    = "../OutputFiles/vascDens_" + to_string(count) + "_" +
                            "1_" + to_string(p) + ".res";
                    nFKilledCells = "../OutputFiles/killedCells_" + to_string(count) + "_" +
                            "1_" + to_string(p) + ".res";
                    nFCycle       = "../OutputFiles/cycle_" + to_string(count) + "_" +
                            "1_" + to_string(p) + ".res";
                    nFHypDens     = "../OutputFiles/hypDens_" + to_string(count) + "_" +
                            "1_" + to_string(p) + ".res";
                    nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(count) + "_" +
                            "1_" + to_string(p) + ".res";
                    nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(count) + "_" +
                            "1_" + to_string(p) + ".res";

                    model(x, y1[p],  nrow, ncol, nlayer, cellSize, inTum, inVes, nFTumDens,
                          nFTumVol, nFVascDens, nFKilledCells, nFCycle, nFHypDens,
                          nFPO2Stat, nFVegfStat);
                    nEv++;

                    cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
                    cout << "---------------------------------------------" << endl;

                    for(int j(0); j < nOut; j++){
                        y1mean[j] += y1[p][j];
                    }
                }
                count++;

                for(int j(0); j < nOut; j++){
                    y1mean[j] /= P;
                }

                fEndTreatTumDens << y0mean[0]  << " " << y1mean[0]  << " ";
                f3MonTumDens     << y0mean[1]  << " " << y1mean[1]  << " ";
                fFinTumVol       << y0mean[2]  << " " << y1mean[2]  << " ";
                fIntTumDens      << y0mean[3]  << " " << y1mean[3]  << " ";
                fKilled50        << y0mean[4]  << " " << y1mean[4]  << " ";
                fKilled80        << y0mean[5]  << " " << y1mean[5]  << " ";
                fKilled90        << y0mean[6]  << " " << y1mean[6]  << " ";
                fKilled95        << y0mean[7]  << " " << y1mean[7]  << " ";
                fTimeTo95        << y0mean[8]  << " " << y1mean[8]  << " ";
                fKilled99        << y0mean[9]  << " " << y1mean[9]  << " ";
                fTimeTo99        << y0mean[10] << " " << y1mean[10] << " ";
                fKilled999       << y0mean[11] << " " << y1mean[11] << " ";
                fRec             << y0mean[12] << " " << y1mean[12] << " ";
                fRecTumDens      << y0mean[13] << " " << y1mean[13] << " ";
                fRecTime         << y0mean[14] << " " << y1mean[14] << " ";

                if(P > 1){
                    for(int j(0); j < nOut; j++){
                        for(int p(0); p < P; p++){
                            y0std[j] += (y0[p][j] - y0mean[j]) * (y0[p][j] - y0mean[j]);
                            y1std[j] += (y1[p][j] - y1mean[j]) * (y1[p][j] - y1mean[j]);
                        }
                        y0std[j] = sqrt(y0std[j] / (P - 1.0));
                        y1std[j] = sqrt(y1std[j] / (P - 1.0));
                    }

                    fEndTreatTumDens << y0std[0]  << " " << y1std[0];
                    f3MonTumDens     << y0std[1]  << " " << y1std[1];
                    fFinTumVol       << y0std[2]  << " " << y1std[2];
                    fIntTumDens      << y0std[3]  << " " << y1std[3];
                    fKilled50        << y0std[4]  << " " << y1std[4];
                    fKilled80        << y0std[5]  << " " << y1std[5];
                    fKilled90        << y0std[6]  << " " << y1std[6];
                    fKilled95        << y0std[7]  << " " << y1std[7];
                    fTimeTo95        << y0std[8]  << " " << y1std[8];
                    fKilled99        << y0std[9]  << " " << y1std[9];
                    fTimeTo99        << y0std[10] << " " << y1std[10];
                    fKilled999       << y0std[11] << " " << y1std[11];
                    fRec             << y0std[12] << " " << y1std[12];
                    fRecTumDens      << y0std[13] << " " << y1std[13];
                    fRecTime         << y0std[14] << " " << y1std[14];
                }

                fEndTreatTumDens << endl;
                f3MonTumDens     << endl;
                fFinTumVol       << endl;
                fIntTumDens      << endl;
                fKilled50        << endl;
                fKilled80        << endl;
                fKilled90        << endl;
                fKilled95        << endl;
                fTimeTo95        << endl;
                fKilled99        << endl;
                fTimeTo99        << endl;
                fKilled999       << endl;
                fRec             << endl;
                fRecTumDens      << endl;
                fRecTime         << endl;

                for (int j(0); j < nMostRelPar; j ++){
                    fCombPar << x[Xi[j]] << " ";
                }
                fCombPar << endl;
            }
        }
    }

    fEndTreatTumDens.close();
    f3MonTumDens.close();
    fFinTumVol.close();
    fIntTumDens.close();
    fKilled50.close();
    fKilled80.close();
    fKilled90.close();
    fKilled95.close();
    fTimeTo95.close();
    fKilled99.close();
    fTimeTo99.close();
    fKilled999.close();
    fRec.close();
    fRecTumDens.close();
    fRecTime.close();
    fCombPar.close();
}


void varParFromFiles(const vector<string> nFPar, const string nFInTissueDim,
                     const string nFInTum, const string nFInVes){
    const int K(39), L(nFPar.size()), nOut(15);
    double x[K], y[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream fEndTreatTumDens, f3MonTumDens, fFinTumVol, fIntTumDens;
    ofstream fKilled50, fKilled80, fKilled90, fKilled95, fKilled99, fKilled999;
    ofstream fTimeTo95, fTimeTo99;
    ofstream fRec, fRecTumDens, fRecTime;

    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                cellSize, inTum, inVes);

    for(int i(0); i < nFPar.size(); i++){
        ifstream fPar(nFPar[i].c_str());
        for(int k(0); k < K; k++){
            fPar >> x[k];
        }
        fPar.close();

        nFTumDens     = "../OutputFiles/tumDens_" + to_string(i) + ".res";
        nFTumVol      = "../OutputFiles/tumVol_" + to_string(i) + ".res";
        nFVascDens    = "../OutputFiles/vascDens_" + to_string(i) + ".res";
        nFKilledCells = "../OutputFiles/killedCells_" + to_string(i) + ".res";
        nFCycle       = "../OutputFiles/cycle_" + to_string(i) + ".res";
        nFHypDens     = "../OutputFiles/hypDens_" + to_string(i) + ".res";
        nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(i) + ".res";
        nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(i) + ".res";

        fEndTreatTumDens.open("../OutputFiles/endTreatTumDens_" + to_string(i) + ".res");
        f3MonTumDens.open("../OutputFiles/3MonTumDens_" + to_string(i) + ".res");
        fFinTumVol.open("../OutputFiles/finTumVol_" + to_string(i) + ".res");
        fIntTumDens.open("../OutputFiles/intTumDens_" + to_string(i) + ".res");
        fKilled50.open("../OutputFiles/killed50_" + to_string(i) + ".res");
        fKilled80.open("../OutputFiles/killed80_" + to_string(i) + ".res");
        fKilled90.open("../OutputFiles/killed90_" + to_string(i) + ".res");
        fKilled95.open("../OutputFiles/killed95_" + to_string(i) + ".res");
        fTimeTo95.open("../OutputFiles/timeTo95_" + to_string(i) + ".res");
        fKilled99.open("../OutputFiles/killed99_" + to_string(i) + ".res");
        fTimeTo99.open("../OutputFiles/timeTo99_" + to_string(i) + ".res");
        fKilled999.open("../OutputFiles/killed99_" + to_string(i) + ".res");
        fRec.open("../OutputFiles/rec_" + to_string(i) + ".res");
        fRecTumDens.open("../OutputFiles/recTumDens_" + to_string(i) + ".res");
        fRecTime.open("../OutputFiles/recTime_" + to_string(i) + ".res");

        model(x, y,  nrow, ncol, nlayer, cellSize, inTum, inVes, nFTumDens,
              nFTumVol, nFVascDens, nFKilledCells, nFCycle, nFHypDens,
              nFPO2Stat, nFVegfStat);

        fEndTreatTumDens << y[0];
        f3MonTumDens     << y[1];
        fFinTumVol       << y[2];
        fIntTumDens      << y[3];
        fKilled50        << y[4];
        fKilled80        << y[5];
        fKilled90        << y[6];
        fKilled95        << y[7];
        fTimeTo95        << y[8];
        fKilled99        << y[9];
        fTimeTo99        << y[10];
        fKilled999       << y[11];
        fRec             << y[12];
        fRecTumDens      << y[13];
        fRecTime         << y[14];

        fEndTreatTumDens.close();
        f3MonTumDens.close();
        fFinTumVol.close();
        fIntTumDens.close();
        fKilled50.close();
        fKilled80.close();
        fKilled90.close();
        fKilled95.close();
        fTimeTo95.close();
        fKilled99.close();
        fTimeTo99.close();
        fKilled999.close();
        fRec.close();
        fRecTumDens.close();
        fRecTime.close();

        cout << i + 1 << " out of " << L << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
    }
}


void varStoch(const int N, const int P, const std::string nFRefParInt,
              const std::string nFInTissueDim, const std::string nFInTum,
              const std::string nFInVes){
    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                cellSize, inTum, inVes);

    const int K(39), nOut(15);
    double h[K], x0[K];
    ifstream fRefParInt(nFRefParInt.c_str());

    for(int k(0); k < K; k++){
        fRefParInt >> x0[k];
        fRefParInt >> h[k];
        h[k] -= x0[k];
    }

    fRefParInt.close();

    double **X, ***Y;

    X = alloc2D(N, K);
    Y = alloc3D(P, N, nOut);

    for(int i(0); i < N; i++){
        for(int k(0); k < K; k++){
            X[i][k] = x0[k] + (double)(rand()) / (double)(RAND_MAX) * h[k];
        }
    }

    ofstream fCombPar("../OutputFiles/combPar.res");

    for(int i(0); i < N; i++){
        for(int k(0); k < K; k++){
            fCombPar << X[i][k] << " ";
        }
        fCombPar << endl;
    }

    int nEv(0);
    const int nEvTot(N * P);

    ofstream fEndTreatTumDens("../OutputFiles/endTreatTumDens.res");
    ofstream f3MonTumDens("../OutputFiles/3MonTumDens.res");
    ofstream fFinTumVol("../OutputFiles/finTumVol.res");
    ofstream fIntTumDens("../OutputFiles/intTumDens.res");
    ofstream fKilled50("../OutputFiles/killed50.res");
    ofstream fKilled80("../OutputFiles/killed80.res");
    ofstream fKilled90("../OutputFiles/killed90.res");
    ofstream fKilled95("../OutputFiles/killed95.res");
    ofstream fTimeTo95("../OutputFiles/timeTo95.res");
    ofstream fKilled99("../OutputFiles/killed99.res");
    ofstream fTimeTo99("../OutputFiles/timeTo99.res");
    ofstream fKilled999("../OutputFiles/killed999.res");
    ofstream fRec("../OutputFiles/rec.res");
    ofstream fRecTumDens("../OutputFiles/recTumDens.res");
    ofstream fRecTime("../OutputFiles/recTime.res");

    for(int i(0); i < N; i++){
        for(int l(0); l < P; l++){
            model(X[i], Y[l][i], nrow, ncol, nlayer, cellSize, inTum, inVes);
            nEv++;
            cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
            cout << "---------------------------------------------" << endl;

            fEndTreatTumDens << Y[l][i][0]  << " ";
            f3MonTumDens     << Y[l][i][1]  << " ";
            fFinTumVol       << Y[l][i][2]  << " ";
            fIntTumDens      << Y[l][i][3]  << " ";
            fKilled50        << Y[l][i][4]  << " ";
            fKilled80        << Y[l][i][5]  << " ";
            fKilled90        << Y[l][i][6]  << " ";
            fKilled95        << Y[l][i][7]  << " ";
            fTimeTo95        << Y[l][i][8]  << " ";
            fKilled99        << Y[l][i][9]  << " ";
            fTimeTo99        << Y[l][i][10] << " ";
            fKilled999       << Y[l][i][11] << " ";
            fRec             << Y[l][i][12] << " ";
            fRecTumDens      << Y[l][i][13] << " ";
            fRecTime         << Y[l][i][14] << " ";
        }

        fEndTreatTumDens << endl;
        f3MonTumDens     << endl;
        fFinTumVol       << endl;
        fIntTumDens      << endl;
        fKilled50        << endl;
        fKilled80        << endl;
        fKilled90        << endl;
        fKilled95        << endl;
        fTimeTo95        << endl;
        fKilled99        << endl;
        fTimeTo99        << endl;
        fKilled999       << endl;
        fRec             << endl;
        fRecTumDens      << endl;
        fRecTime         << endl;
    }

    fEndTreatTumDens.close();
    f3MonTumDens.close();
    fFinTumVol.close();
    fIntTumDens.close();
    fKilled50.close();
    fKilled80.close();
    fKilled90.close();
    fKilled95.close();
    fTimeTo95.close();
    fKilled99.close();
    fTimeTo99.close();
    fKilled999.close();
    fRec.close();
    fRecTumDens.close();
    fRecTime.close();

    free2D(X, N);
    free3D(Y, P, N);
}

