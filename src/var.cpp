#include "sensAn.hpp"

using namespace std;

void var1ParRange(const int kp, const int L, const string nRefParInt,
                  const string nFInTissueDim, const string nFInTum,
                  const string nFInVes){
    const int K(38);
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
            const string nFInTissueDim, const string nFInTum, const string nFInVes){
    const int K(38), nOut(15), p(6);
    int nLeastRelPar, nMostRelPar, nVarPar;

    ifstream fMostRelPar(nFMostRelPar.c_str());

    fMostRelPar >> nMostRelPar;

    int Xi[nMostRelPar];
    double h, x[K], X[nMostRelPar][p];

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

    int nEv(0), nEvTot(2 * pow(p, nMostRelPar));
    double err[nOut], errRel[nOut], maxy[nOut], y0[nOut], y1[nOut];
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

    for(int i1(0); i1 < p; i1++){
        x[Xi[0]] = X[0][i1];
        for(int i2(0); i2 < p; i2++){
            x[Xi[1]] = X[1][i2];
            for(int i3(0); i3 < p; i3++){
                x[Xi[2]] = X[2][i3];
                for(int k(0); k < nVarPar; k++){
                    x[XVari[k]] = XVar[k][0];
                }

                nFTumDens     = "../OutputFiles/tumDens_" + to_string(nEv) + ".res";
                nFTumVol      = "../OutputFiles/tumVol_" + to_string(nEv) + ".res";
                nFVascDens    = "../OutputFiles/vascDens_" + to_string(nEv) + ".res";
                nFKilledCells = "../OutputFiles/killedCells_" + to_string(nEv) + ".res";
                nFCycle       = "../OutputFiles/cycle_" + to_string(nEv) + ".res";
                nFHypDens     = "../OutputFiles/hypDens_" + to_string(nEv) + ".res";
                nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(nEv) + ".res";
                nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(nEv) + ".res";

                model(x, y0,  nrow, ncol, nlayer, cellSize, inTum, inVes, nFTumDens,
                      nFTumVol, nFVascDens, nFKilledCells, nFCycle, nFHypDens,
                      nFPO2Stat, nFVegfStat);
                nEv++;

                cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
                cout << "---------------------------------------------" << endl;

                for(int k(0); k < nVarPar; k++){
                    x[XVari[k]] = XVar[k][1];
                }

                nFTumDens     = "../OutputFiles/tumDens_" + to_string(nEv) + ".res";
                nFTumVol      = "../OutputFiles/tumVol_" + to_string(nEv) + ".res";
                nFVascDens    = "../OutputFiles/vascDens_" + to_string(nEv) + ".res";
                nFKilledCells = "../OutputFiles/killedCells_" + to_string(nEv) + ".res";
                nFCycle       = "../OutputFiles/cycle_" + to_string(nEv) + ".res";
                nFHypDens     = "../OutputFiles/hypDens_" + to_string(nEv) + ".res";
                nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(nEv) + ".res";
                nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(nEv) + ".res";

                model(x, y1,  nrow, ncol, nlayer, cellSize, inTum, inVes, nFTumDens,
                      nFTumVol, nFVascDens, nFKilledCells, nFCycle, nFHypDens,
                      nFPO2Stat, nFVegfStat);
                nEv++;

                cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
                cout << "---------------------------------------------" << endl;

                for(int j(0); j < nOut; j++){
                    err[j] = fabs(y0[j] - y1[j]);
                    maxy[j] = max(y0[j], y1[j]);
                    if(maxy[j]){
                        errRel[j] = err[j] / maxy[j];
                    }
                    else{
                        errRel[j] = 0.0;
                    }
                }

                fEndTreatTumDens << y0[0]  << " " << y1[0] << endl;
                f3MonTumDens     << y0[1]  << " " << y1[1] << endl;
                fFinTumVol       << y0[2]  << " " << y1[2] << endl;
                fIntTumDens      << y0[3]  << " " << y1[3] << endl;
                fKilled50        << y0[4]  << " " << y1[4] << endl;
                fKilled80        << y0[5]  << " " << y1[5] << endl;
                fKilled90        << y0[6]  << " " << y1[6] << endl;
                fKilled95        << y0[7]  << " " << y1[7] << endl;
                fTimeTo95        << y0[8]  << " " << y1[8] << endl;
                fKilled99        << y0[9]  << " " << y1[9] << endl;
                fTimeTo99        << y0[10] << " " << y1[10] << endl;
                fKilled999       << y0[11] << " " << y1[11] << endl;
                fRec             << y0[12] << " " << y1[12] << endl;
                fRecTumDens      << y0[13] << " " << y1[13] << endl;
                fRecTime         << y0[14] << " " << y1[14] << endl;

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
    const int K(38), L(nFPar.size()), nOut(15);
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
