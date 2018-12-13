#include "sensAn.hpp"

using namespace std;

void var1ParRange(const int kp, const int L, const string nRefParInt){
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

    const int nOut(8);
    const double delta(h[kp] / (L - 1));
    double y[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream fEndTreatTumDens, f3MonTumDens, fRecTumDens;
    ofstream fFinTumVol, fIntTumDens, fTimeTo95;
    ofstream fTimeTo99, fRecTime;

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
        fRecTumDens.open("../OutputFiles/recTumDens_" + to_string(i) + ".res");
        fFinTumVol.open("../OutputFiles/finTumVol_" + to_string(i) + ".res");
        fIntTumDens.open("../OutputFiles/intTumDens_" + to_string(i) + ".res");
        fTimeTo95.open("../OutputFiles/timeTo95_" + to_string(i) + ".res");
        fTimeTo99.open("../OutputFiles/timeTo99_" + to_string(i) + ".res");
        fRecTime.open("../OutputFiles/recTime_" + to_string(i) + ".res");

        model(x, y, nFTumDens, nFTumVol, nFVascDens, nFKilledCells,
              nFCycle, nFHypDens, nFPO2Stat, nFVegfStat);

        fEndTreatTumDens << y[0];
        f3MonTumDens     << y[1];
        fRecTumDens      << y[2];
        fFinTumVol       << y[3];
        fIntTumDens      << y[4];
        fTimeTo95        << y[5];
        fTimeTo99        << y[6];
        fRecTime         << y[7];

        fEndTreatTumDens.close();
        f3MonTumDens.close();
        fRecTumDens.close();
        fFinTumVol.close();
        fIntTumDens.close();
        fTimeTo95.close();
        fTimeTo99.close();
        fRecTime.close();

        cout << i + 1 << " out of " << L << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
        x[kp] += delta;
    }
}


void varErr(const string nFVarPar, const string nFMostRelPar, const string nFLeastPar){
    const int K(38), nOut(8), p(6);
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

    int nEv(0), nEvTot(2 * pow(p, nMostRelPar));
    double err[nOut], errRel[nOut], maxy[nOut], y0[nOut], y1[nOut];
    ofstream fErrEndTreatTumDens("../OutputFiles/errEndTreatTumDens.res");
    ofstream fErr3MonTumDens("../OutputFiles/err3MonTumDens.res");
    ofstream fErrRecTumDens("../OutputFiles/errRecTumDens.res");
    ofstream fErrFinTumVol("../OutputFiles/errFinTumVol.res");
    ofstream fErrIntTumDens("../OutputFiles/errIntTumDens.res");
    ofstream fErrTimeTo95("../OutputFiles/errTimeTo95.res");
    ofstream fErrTimeTo99("../OutputFiles/errTimeTo99.res");
    ofstream fErrRecTime("../OutputFiles/errRecTime.res");

    for(int i1(0); i1 < p; i1++){
        x[Xi[0]] = X[0][i1];
        for(int i2(0); i2 < p; i2++){
            x[Xi[1]] = X[1][i2];
            for(int i3(0); i3 < p; i3++){
                x[Xi[2]] = X[2][i3];
                for(int i4(0); i4 < p; i4++){
                    x[Xi[3]] = X[3][i4];
                    for(int k(0); k < nVarPar; k++){
                        x[XVari[k]] = XVar[k][0];
                    }
                    model(x, y0);
                    nEv++;

                    cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
                    cout << "---------------------------------------------" << endl;

                    for(int k(0); k < nVarPar; k++){
                        x[XVari[k]] = XVar[k][1];
                    }
                    model(x, y1);
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

                    fErrEndTreatTumDens << err[0] << " " << errRel[0] << " ";
                    fErr3MonTumDens     << err[1] << " " << errRel[1] << " ";
                    fErrRecTumDens      << err[2] << " " << errRel[2] << " ";
                    fErrFinTumVol       << err[3] << " " << errRel[3] << " ";
                    fErrIntTumDens      << err[4] << " " << errRel[4] << " ";
                    fErrTimeTo95        << err[5] << " " << errRel[5] << " ";
                    fErrTimeTo99        << err[6] << " " << errRel[6] << " ";
                    fErrRecTime         << err[7] << " " << errRel[7] << " ";

                    for (int j(0); j < nMostRelPar; j ++){
                        fErrEndTreatTumDens << x[Xi[j]] << " ";
                        fErr3MonTumDens     << x[Xi[j]] << " ";
                        fErrRecTumDens      << x[Xi[j]] << " ";
                        fErrFinTumVol       << x[Xi[j]] << " ";
                        fErrIntTumDens      << x[Xi[j]] << " ";
                        fErrTimeTo95        << x[Xi[j]] << " ";
                        fErrTimeTo99        << x[Xi[j]] << " ";
                        fErrRecTime         << x[Xi[j]] << " ";
                    }
                    fErrEndTreatTumDens << endl;
                    fErr3MonTumDens     << endl;
                    fErrRecTumDens      << endl;
                    fErrFinTumVol       << endl;
                    fErrIntTumDens      << endl;
                    fErrTimeTo95        << endl;
                    fErrTimeTo99        << endl;
                    fErrRecTime         << endl;

                }
            }
        }
    }

    fErrEndTreatTumDens.close();
    fErr3MonTumDens.close();
    fErrRecTumDens.close();
    fErrFinTumVol.close();
    fErrIntTumDens.close();
    fErrTimeTo95.close();
    fErrTimeTo99.close();
    fErrRecTime.close();
}

void varParFromFiles(const vector<string> nFPar){
    const int K(38), L(nFPar.size()), nOut(8);
    double x[K], y[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream fEndTreatTumDens, f3MonTumDens, fRecTumDens;
    ofstream fFinTumVol, fIntTumDens, fTimeTo95;
    ofstream fTimeTo99, fRecTime;

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
        fRecTumDens.open("../OutputFiles/recTumDens_" + to_string(i) + ".res");
        fFinTumVol.open("../OutputFiles/finTumVol_" + to_string(i) + ".res");
        fIntTumDens.open("../OutputFiles/intTumDens_" + to_string(i) + ".res");
        fTimeTo95.open("../OutputFiles/timeTo95_" + to_string(i) + ".res");
        fTimeTo99.open("../OutputFiles/timeTo99_" + to_string(i) + ".res");
        fRecTime.open("../OutputFiles/recTime_" + to_string(i) + ".res");

        model(x, y, nFTumDens, nFTumVol, nFVascDens, nFKilledCells,
              nFCycle, nFHypDens, nFPO2Stat, nFVegfStat);

        fEndTreatTumDens << y[0];
        f3MonTumDens     << y[1];
        fRecTumDens      << y[2];
        fFinTumVol       << y[3];
        fIntTumDens      << y[4];
        fTimeTo95        << y[5];
        fTimeTo99        << y[6];
        fRecTime         << y[7];

        fEndTreatTumDens.close();
        f3MonTumDens.close();
        fRecTumDens.close();
        fFinTumVol.close();
        fIntTumDens.close();
        fTimeTo95.close();
        fTimeTo99.close();
        fRecTime.close();

        cout << i + 1 << " out of " << L << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
    }
}
