#include "var.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * This function studies the impact of incrementing the value of one parameter
 * on the scalar outputs of the model. Regular increments are considered. Mean
 * values of ranges are used for the not varying parameters.
 *
 * Inputs:
 *  - kp: parameter of the model to be incremented,
 *  - L: number of regular increments,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration,
 *  - nFInPO2: name of the file containing the initial pO2 values.
------------------------------------------------------------------------------*/

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

    const int nOut(20);
    const double delta(h[kp] / (L - 1));
    double y[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream f8wTumDens, f12wTumDens, f8wTumVol, f12wTumVol;
    ofstream f8wIntTumDens, f12wIntTumDens, f8wIntTumVol, f12wIntTumVol;
    ofstream fKilled50, fKilled80, fKilled90, fKilled95, fKilled99, fKilled999;
    ofstream fTimeTo95, fTimeTo99;
    ofstream fRec, fRecTumDens, fRecTumVol, fRecTime;

    for(int i(0); i < L; i++){
        nFTumDens     = "../OutputFiles/tumDens_" + to_string(i) + ".res";
        nFTumVol      = "../OutputFiles/tumVol_" + to_string(i) + ".res";
        nFVascDens    = "../OutputFiles/vascDens_" + to_string(i) + ".res";
        nFKilledCells = "../OutputFiles/killedCells_" + to_string(i) + ".res";
        nFDeadDens    = "../OutputFiles/deadDens_" + to_string(i) + ".res";
        nFCycle       = "../OutputFiles/cycle_" + to_string(i) + ".res";
        nFHypDens     = "../OutputFiles/hypDens_" + to_string(i) + ".res";
        nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(i) + ".res";
        nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(i) + ".res";

        f8wTumDens.open("../OutputFiles/8wTumDens_" + to_string(i) +
                              ".res");
        f12wTumDens.open("../OutputFiles/12wTumDens_" + to_string(i) +
                          ".res");
        f8wTumVol.open("../OutputFiles/8wTumVol_" + to_string(i) +
                              ".res");
        f12wTumVol.open("../OutputFiles/12wTumVol_" + to_string(i) +
                          ".res");
        f8wIntTumDens.open("../OutputFiles/8wIntTumDens_" + to_string(i) +
                              ".res");
        f12wIntTumDens.open("../OutputFiles/12wIntTumDens_" + to_string(i) +
                          ".res");
        f8wIntTumVol.open("../OutputFiles/8wIntTumVol_" + to_string(i) +
                              ".res");
        f12wIntTumVol.open("../OutputFiles/12wIntTumVol_" + to_string(i) +
                          ".res");
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
        fRecTumVol.open("../OutputFiles/recTumVol_" + to_string(i) + ".res");
        fRecTime.open("../OutputFiles/recTime_" + to_string(i) + ".res");

        model(x, y, nrow, ncol, nlayer, cellSize, inTum, inVes, nFTumDens,
              nFTumVol, nFVascDens, nFKilledCells, nFDeadDens, nFCycle,
              nFHypDens, nFPO2Stat, nFVegfStat);

        f8wTumDens     << y[0];
        f12wTumDens    << y[1];
        f8wTumVol      << y[2];
        f12wTumVol     << y[3];
        f8wIntTumDens  << y[4];
        f12wIntTumDens << y[5];
        f8wIntTumVol   << y[6];
        f12wIntTumVol  << y[7];
        fKilled50      << y[8];
        fKilled80      << y[9];
        fKilled90      << y[10];
        fKilled95      << y[11];
        fTimeTo95      << y[12];
        fKilled99      << y[13];
        fTimeTo99      << y[14];
        fKilled999     << y[15];
        fRec           << y[16];
        fRecTumDens    << y[17];
        fRecTumVol     << y[18];
        fRecTime       << y[19];

        f8wTumDens.close();
        f12wTumDens.close();
        f8wTumVol.close();
        f12wTumVol.close();
        f8wIntTumDens.close();
        f12wIntTumDens.close();
        f8wIntTumVol.close();
        f12wIntTumVol.close();
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
        fRecTumVol.close();
        fRecTime.close(););
        fRecTime.close();

        cout << i + 1 << " out of " << L << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
        x[kp] += delta;
    }
}


/*------------------------------------------------------------------------------
 * This function studies the impact of the initial tumour and vascular
 * architectures. Mean values of ranges are used for the parameters.
 *
 * Inputs:
 *  - P: number of repetitions of each combination,
 *  - nFArt: name of the file containing the values of the tumour and
 *  vascular densities and sigmas,
 *  - nFInTissueDim: name of the file containing the dimensions of an artificial
 *  tissue,
 *  - nFRefParMean: name of the file containing the mean values of the
 *  parameters.
------------------------------------------------------------------------------*/

void varArtTissue(const int P, const string nFArt,
                  const string nFInTissueDim, const string nFRefParMean){
    const int K(39), nOut(20);
    int nrow, ncol, nlayer;
    double cellSize;
    readInFiles(nFInTissueDim, nrow, ncol, nlayer, cellSize);

    int nTumDens, nSigmaTum, nVascDens, nSigmaVasc;
    ifstream fArt(nFArt.c_str());
    fArt >> nTumDens >> nSigmaTum >> nVascDens >> nSigmaVasc;

    double tumDens[nTumDens], sigmaTum[nSigmaTum];
    double vascDens[nVascDens], sigmaVasc[nSigmaVasc];

    for(int i(0); i < nTumDens; i++){
        fArt >> tumDens[i];
    }
    for(int i(0); i < nSigmaTum; i++){
        fArt >> sigmaTum[i];
    }
    for(int i(0); i < nVascDens; i++){
        fArt >> vascDens[i];
    }
    for(int i(0); i < nSigmaVasc; i++){
        fArt >> sigmaVasc[i];
    }
    fArt.close();

    double x[K];
    ifstream fRefParMean(nFRefParMean.c_str());

    for(int k(0); k < K; k++){
        fRefParMean >> x[k];
    }
    fRefParMean.close();

    const int NN(nTumDens * nSigmaTum * nVascDens * nSigmaVasc);
    const int nEvTot(NN * P), nrowNcolNlayer(nrow * ncol *nlayer);
    int i(0), nEv(0);
    double y[P][nOut], ymean[nOut], ystd[nOut];
    vector<bool> inTum(nrowNcolNlayer), inVes(nrowNcolNlayer);
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream f8wTumDens("../OutputFiles/8wTumDens.res");
    ofstream f12wTumDens("../OutputFiles/12wTumDens.res");
    ofstream f8wTumVol("../OutputFiles/8wTumVol.res");
    ofstream f12wTumVol("../OutputFiles/12wTumVol.res");
    ofstream f8wIntTumDens("../OutputFiles/8wIntTumDens.res");
    ofstream f12wIntTumDens("../OutputFiles/12wIntTumDens.res");
    ofstream f8wIntTumVol("../OutputFiles/8wIntTumVol.res");
    ofstream f12wIntTumVol("../OutputFiles/12wIntTumVol.res");
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
    ofstream fRecTumVol("../OutputFiles/recTumVol.res");
    ofstream fRecTime("../OutputFiles/recTime.res");
    ofstream fComb("../OutputFiles/comb.res");

    for(int iTumDens(0); iTumDens < nTumDens; iTumDens++){
        for(int iSigmaTum(0); iSigmaTum < nSigmaTum; iSigmaTum++){
            for(int iVascDens(0); iVascDens < nVascDens; iVascDens++){
                for(int iSigmaVasc(0); iSigmaVasc < nSigmaVasc; iSigmaVasc++){
                    for(int j(0); j < nOut; j++){
                        ymean[j] = 0.0;
                        ystd[j] = 0.0;
                    }
                    for(int p(0); p < P; p++){
                        nFTumDens     = "../OutputFiles/tumDens_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        nFTumVol      = "../OutputFiles/tumVol_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        nFVascDens    = "../OutputFiles/vascDens_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        nFKilledCells = "../OutputFiles/killedCells_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        nFDeadDens    = "../OutputFiles/deadDens_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        nFCycle       = "../OutputFiles/cycle_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        nFHypDens     = "../OutputFiles/hypDens_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        nFPO2Stat     = "../OutputFiles/pO2Stat_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        nFVegfStat    = "../OutputFiles/vegfStat_" +
                                to_string(i) + "_" + to_string(p) + ".res";
                        cout << tumDens[iTumDens]      << " " <<
                                sigmaTum[iSigmaTum]   << " " <<
                                vascDens[iVascDens]   << " " <<
                                sigmaVasc[iSigmaVasc] << endl;
                        createInFiles(nrow, ncol, nlayer, tumDens[iTumDens],
                                      sigmaTum[iSigmaTum], vascDens[iVascDens],
                                      sigmaVasc[iSigmaVasc], inTum, inVes);
                        model(x, y[p], nrow, ncol, nlayer, cellSize, inTum,
                              inVes, nFTumDens, nFTumVol, nFVascDens,
                              nFKilledCells, nFDeadDens, nFCycle, nFHypDens,
                              nFPO2Stat, nFVegfStat);
                        nEv++;

                        cout << nEv << " out of " << nEvTot <<
                                " evaluations of the model" << endl;
                        cout << "-------------------------------------" << endl;

                        for(int j(0); j < nOut; j++){
                            ymean[j] += y[p][j];
                        }
                    }

                    for(int j(0); j < nOut; j++){
                        ymean[j] /= P;
                    }

                    f8wTumDens     << ymean[0]   << " ";
                    f12wTumDens    << ymean[1]   << " ";
                    f8wTumVol      << ymean[2]   << " ";
                    f12wTumVol     << ymean[3]   << " ";
                    f8wIntTumDens  << ymean[4]   << " ";
                    f12wIntTumDens << ymean[5]   << " ";
                    f8wIntTumVol   << ymean[6]   << " ";
                    f12wIntTumVol  << ymean[7]   << " ";
                    fKilled50      << ymean[8]   << " ";
                    fKilled80      << ymean[9]   << " ";
                    fKilled90      << ymean[10]  << " ";
                    fKilled95      << ymean[11]  << " ";
                    fTimeTo95      << ymean[12]  << " ";
                    fKilled99      << ymean[13]  << " ";
                    fTimeTo99      << ymean[14]  << " ";
                    fKilled999     << ymean[15]  << " ";
                    fRec           << ymean[16]  << " ";
                    fRecTumDens    << ymean[17]  << " ";
                    fRecTumVol     << ymean[18]  << " ";
                    fRecTime       << ymean[19]  << " ";

                    if(P > 1){
                        for(int j(0); j < nOut; j++){
                            for(int p(0); p < P; p++){
                                ystd[j] += (y[p][j] - ymean[j]) * (y[p][j] -
                                                                   ymean[j]);
                            }
                            ystd[j] = sqrt(ystd[j] / (P - 1.0));
                        }

                        f8wTumDens     << ystd[0]   << " ";
                        f12wTumDens    << ystd[1]   << " ";
                        f8wTumVol      << ystd[2]   << " ";
                        f12wTumVol     << ystd[3]   << " ";
                        f8wIntTumDens  << ystd[4]   << " ";
                        f12wIntTumDens << ystd[5]   << " ";
                        f8wIntTumVol   << ystd[6]   << " ";
                        f12wIntTumVol  << ystd[7]   << " ";
                        fKilled50      << ystd[8]   << " ";
                        fKilled80      << ystd[9]   << " ";
                        fKilled90      << ystd[10]  << " ";
                        fKilled95      << ystd[11]  << " ";
                        fTimeTo95      << ystd[12]  << " ";
                        fKilled99      << ystd[13]  << " ";
                        fTimeTo99      << ystd[14]  << " ";
                        fKilled999     << ystd[15]  << " ";
                        fRec           << ystd[16]  << " ";
                        fRecTumDens    << ystd[17]  << " ";
                        fRecTumVol     << ystd[18]  << " ";
                        fRecTime       << ystd[19]  << " ";
                    }

                    f8wTumDens     << endl;
                    f12wTumDens    << endl;
                    f8wTumVol      << endl;
                    f12wTumVol     << endl;
                    f8wIntTumDens  << endl;
                    f12wIntTumDens << endl;
                    f8wIntTumVol   << endl;
                    f12wIntTumVol  << endl;
                    fKilled50      << endl;
                    fKilled80      << endl;
                    fKilled90      << endl;
                    fKilled95      << endl;
                    fTimeTo95      << endl;
                    fKilled99      << endl;
                    fTimeTo99      << endl;
                    fKilled999     << endl;
                    fRec           << endl;
                    fRecTumDens    << endl;
                    fRecTumVol     << endl;
                    fRecTime       << endl;

                    fComb << tumDens[iTumDens]     << " " <<
                             sigmaTum[iSigmaTum]   << " " <<
                             vascDens[iVascDens]   << " " <<
                             sigmaVasc[iSigmaVasc] << endl;
                    i++;
                }
            }
        }
    }

    f8wTumDens.close();
    f12wTumDens.close();
    f8wTumVol.close();
    f12wTumVol.close();
    f8wIntTumDens.close();
    f12wIntTumDens.close();
    f8wIntTumVol.close();
    f12wIntTumVol.close();
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
    fRecTumVol.close();
    fRecTime.close(););
    fRecTime.close();
    fComb.close();
}


/*------------------------------------------------------------------------------
 * This function studies the impact of the initial tumour and vascular
 * densities, supposing a random uniform distribution. Mean values of ranges are
 * used for the parameters.
 *
 * Inputs:
 *  - N: number of possible values of each density (N x N combinations, in
 *  total),
 *  - P: number of repetitions of each combination,
 *  - nFDensInt: name of the file containing the intervals of the tumour and
 *  vascular densities,
 *  - nFInTissueDim: name of the file containing the dimensions of an artificial
 *  tissue,
 *  - nFRefParMean: name of the file containing the mean values of the
 *  parameters.
------------------------------------------------------------------------------*/

void varArtTissue(const int N, const int P, const string nFDensInt,
                  const string nFInTissueDim, const string nFRefParMean){
    const int K(39), nOut(20);
    int nrow, ncol, nlayer;
    double cellSize;
    readInFiles(nFInTissueDim, nrow, ncol, nlayer, cellSize);

    double tumDensMin, tumDensMax, vascDensMin, vascDensMax;
    ifstream fDensInt(nFDensInt.c_str());

    fDensInt >> tumDensMin >> tumDensMax;
    fDensInt >> vascDensMin >> vascDensMax;
    fDensInt.close();

    double hTumDens, hVascDens;
    hTumDens = (tumDensMax - tumDensMin) / (N - 1);
    hVascDens = (vascDensMax - vascDensMin) / (N - 1);

    double x[K];
    ifstream fRefParMean(nFRefParMean.c_str());

    for(int k(0); k < K; k++){
        fRefParMean >> x[k];
    }
    fRefParMean.close();

    const int nEvTot(N * N * P), nrowNcolNlayer(nrow * ncol *nlayer);
    int i(0), nEv(0);
    double tumDens(tumDensMin), vascDens(vascDensMin);
    double y[P][nOut], ymean[nOut], ystd[nOut];
    vector<bool> inTum(nrowNcolNlayer), inVes(nrowNcolNlayer);
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream f8wTumDens("../OutputFiles/8wTumDens.res");
    ofstream f12wTumDens("../OutputFiles/12wTumDens.res");
    ofstream f8wTumVol("../OutputFiles/8wTumVol.res");
    ofstream f12wTumVol("../OutputFiles/12wTumVol.res");
    ofstream f8wIntTumDens("../OutputFiles/8wIntTumDens.res");
    ofstream f12wIntTumDens("../OutputFiles/12wIntTumDens.res");
    ofstream f8wIntTumVol("../OutputFiles/8wIntTumVol.res");
    ofstream f12wIntTumVol("../OutputFiles/12wIntTumVol.res");
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
    ofstream fRecTumVol("../OutputFiles/recTumVol.res");
    ofstream fRecTime("../OutputFiles/recTime.res");
    ofstream fCombDens("../OutputFiles/combDens.res");

    for(int iTum(0); iTum < N; iTum++){
        vascDens = vascDensMin;
        for(int iVasc(0); iVasc < N; iVasc++){
            for(int j(0); j < nOut; j++){
                ymean[j] = 0.0;
                ystd[j] = 0.0;
            }

            for(int p(0); p < P; p++){
                nFTumDens     = "../OutputFiles/tumDens_" + to_string(i) + "_" +
                        to_string(p) + ".res";
                nFTumVol      = "../OutputFiles/tumVol_" + to_string(i) + "_" +
                        to_string(p) + ".res";
                nFVascDens    = "../OutputFiles/vascDens_" + to_string(i) +
                        "_" + to_string(p) + ".res";
                nFKilledCells = "../OutputFiles/killedCells_" + to_string(i) +
                        "_" + to_string(p) + ".res";
                nFDeadDens    = "../OutputFiles/deadDens_" + to_string(i) +
                        "_" + to_string(p) + ".res";
                nFCycle       = "../OutputFiles/cycle_" + to_string(i) + "_" +
                        to_string(p) + ".res";
                nFHypDens     = "../OutputFiles/hypDens_" + to_string(i) + "_" +
                        to_string(p) + ".res";
                nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(i) + "_" +
                        to_string(p) + ".res";
                nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(i) + "_"
                        + to_string(p) + ".res";
                createInFiles(nrow, ncol, nlayer, tumDens, vascDens, inTum,
                              inVes);
                model(x, y[p], nrow, ncol, nlayer, cellSize, inTum, inVes,
                      nFTumDens, nFTumVol, nFVascDens, nFKilledCells,
                      nFDeadDens, nFCycle, nFHypDens, nFPO2Stat, nFVegfStat);
                nEv++;

                cout << nEv << " out of " << nEvTot <<
                        " evaluations of the model" << endl;
                cout << "---------------------------------------------" << endl;


                for(int j(0); j < nOut; j++){
                    ymean[j] += y[p][j];
                }
            }

            for(int j(0); j < nOut; j++){
                ymean[j] /= P;
            }

            f8wTumDens     << ymean[0]   << " ";
            f12wTumDens    << ymean[1]   << " ";
            f8wTumVol      << ymean[2]   << " ";
            f12wTumVol     << ymean[3]   << " ";
            f8wIntTumDens  << ymean[4]   << " ";
            f12wIntTumDens << ymean[5]   << " ";
            f8wIntTumVol   << ymean[6]   << " ";
            f12wIntTumVol  << ymean[7]   << " ";
            fKilled50      << ymean[8]   << " ";
            fKilled80      << ymean[9]   << " ";
            fKilled90      << ymean[10]  << " ";
            fKilled95      << ymean[11]  << " ";
            fTimeTo95      << ymean[12]  << " ";
            fKilled99      << ymean[13]  << " ";
            fTimeTo99      << ymean[14]  << " ";
            fKilled999     << ymean[15]  << " ";
            fRec           << ymean[16]  << " ";
            fRecTumDens    << ymean[17]  << " ";
            fRecTumVol     << ymean[18]  << " ";
            fRecTime       << ymean[19]  << " ";

            if(P > 1){
                for(int j(0); j < nOut; j++){
                    for(int p(0); p < P; p++){
                        ystd[j] += (y[p][j] - ymean[j]) * (y[p][j] - ymean[j]);
                    }
                    ystd[j] = sqrt(ystd[j] / (P - 1.0));
                }

                f8wTumDens     << ystd[0]   << " ";
                f12wTumDens    << ystd[1]   << " ";
                f8wTumVol      << ystd[2]   << " ";
                f12wTumVol     << ystd[3]   << " ";
                f8wIntTumDens  << ystd[4]   << " ";
                f12wIntTumDens << ystd[5]   << " ";
                f8wIntTumVol   << ystd[6]   << " ";
                f12wIntTumVol  << ystd[7]   << " ";
                fKilled50      << ystd[8]   << " ";
                fKilled80      << ystd[9]   << " ";
                fKilled90      << ystd[10]  << " ";
                fKilled95      << ystd[11]  << " ";
                fTimeTo95      << ystd[12]  << " ";
                fKilled99      << ystd[13]  << " ";
                fTimeTo99      << ystd[14]  << " ";
                fKilled999     << ystd[15]  << " ";
                fRec           << ystd[16]  << " ";
                fRecTumDens    << ystd[17]  << " ";
                fRecTumVol     << ystd[18]  << " ";
                fRecTime       << ystd[19]  << " ";
            }

            f8wTumDens     << endl;
            f12wTumDens    << endl;
            f8wTumVol      << endl;
            f12wTumVol     << endl;
            f8wIntTumDens  << endl;
            f12wIntTumDens << endl;
            f8wIntTumVol   << endl;
            f12wIntTumVol  << endl;
            fKilled50      << endl;
            fKilled80      << endl;
            fKilled90      << endl;
            fKilled95      << endl;
            fTimeTo95      << endl;
            fKilled99      << endl;
            fTimeTo99      << endl;
            fKilled999     << endl;
            fRec           << endl;
            fRecTumDens    << endl;
            fRecTumVol     << endl;
            fRecTime       << endl;

            fCombDens << tumDens << " " << vascDens << endl;
            vascDens += hVascDens;
            i++;
        }
        tumDens += hTumDens;
    }

    f8wTumDens.close();
    f12wTumDens.close();
    f8wTumVol.close();
    f12wTumVol.close();
    f8wIntTumDens.close();
    f12wIntTumDens.close();
    f8wIntTumVol.close();
    f12wIntTumVol.close();
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
    fRecTumVol.close();
    fRecTime.close(););
    fRecTime.close();
    fCombDens.close();
}


/*------------------------------------------------------------------------------
 * This function studies the impact of the initial tumour and vascular
 * architecture, supposing a random uniform distribution of endothelial cells.
 * Mean values of ranges are used for the parameters.
 *
 * Inputs:
 *  - N: number of possible values of each density (N x N combinations, in
 *  total),
 *  - P: number of repetitions of each combination,
 *  - nFDensInt: name of the file containing the intervals of the tumour and
 *  vascular densities,
 *  - nFSigmaTumInt: name of the file containing the intervals of the tumour
 *  sigma,
 *  - nFInTissueDim: name of the file containing the dimensions of an artificial
 *  tissue,
 *  - nFRefParMean: name of the file containing the mean values of the
 *  parameters.
------------------------------------------------------------------------------*/

void varArtTissue(const int N, const int P, const string nFDensInt,
                  const string nFSigmaTumInt, const string nFInTissueDim,
                  const string nFRefParMean){
    const int K(39), nOut(20);
    int nrow, ncol, nlayer;
    double cellSize;
    readInFiles(nFInTissueDim, nrow, ncol, nlayer, cellSize);

    double tumDensMin, tumDensMax, vascDensMin, vascDensMax;
    ifstream fDensInt(nFDensInt.c_str());

    fDensInt >> tumDensMin >> tumDensMax;
    fDensInt >> vascDensMin >> vascDensMax;
    fDensInt.close();

    double sigmaTumMin, sigmaTumMax;
    ifstream fSigmaTumInt(nFSigmaTumInt.c_str());

    fSigmaTumInt >> sigmaTumMin >> sigmaTumMax;
    fSigmaTumInt.close();

    double hTumDens, hSigmaTum, hVascDens;
    hTumDens = (tumDensMax - tumDensMin) / (N - 1);
    hSigmaTum = (sigmaTumMax - sigmaTumMin) / (N - 1);
    hVascDens = (vascDensMax - vascDensMin) / (N - 1);

    double x[K];
    ifstream fRefParMean(nFRefParMean.c_str());

    for(int k(0); k < K; k++){
        fRefParMean >> x[k];
    }
    fRefParMean.close();

    const int nEvTot(N * N * N * P), nrowNcolNlayer(nrow * ncol *nlayer);
    int i(0), nEv(0);
    double tumDens(tumDensMin), sigmaTum(sigmaTumMin), vascDens(vascDensMin);
    double y[P][nOut], ymean[nOut], ystd[nOut];
    vector<bool> inTum(nrowNcolNlayer), inVes(nrowNcolNlayer);
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream f8wTumDens("../OutputFiles/8wTumDens.res");
    ofstream f12wTumDens("../OutputFiles/12wTumDens.res");
    ofstream f8wTumVol("../OutputFiles/8wTumVol.res");
    ofstream f12wTumVol("../OutputFiles/12wTumVol.res");
    ofstream f8wIntTumDens("../OutputFiles/8wIntTumDens.res");
    ofstream f12wIntTumDens("../OutputFiles/12wIntTumDens.res");
    ofstream f8wIntTumVol("../OutputFiles/8wIntTumVol.res");
    ofstream f12wIntTumVol("../OutputFiles/12wIntTumVol.res");
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
    ofstream fRecTumVol("../OutputFiles/recTumVol.res");
    ofstream fRecTime("../OutputFiles/recTime.res");
    ofstream fCombDens("../OutputFiles/combDens.res");

    for(int iTumDens(0); iTumDens < N; iTumDens++){
        sigmaTum = sigmaTumMin;
        for(int iSigmaTum(0); iSigmaTum < N; iSigmaTum++){;
            vascDens = vascDensMin;
            for(int iVasc(0); iVasc < N; iVasc++){
                for(int j(0); j < nOut; j++){
                    ymean[j] = 0.0;
                    ystd[j] = 0.0;
                }

                for(int p(0); p < P; p++){
                    nFTumDens     = "../OutputFiles/tumDens_" + to_string(i) +
                            "_" + to_string(p) + ".res";
                    nFTumVol      = "../OutputFiles/tumVol_" + to_string(i) +
                            "_" + to_string(p) + ".res";
                    nFVascDens    = "../OutputFiles/vascDens_" + to_string(i) +
                            "_" + to_string(p) + ".res";
                    nFKilledCells = "../OutputFiles/killedCells_" +
                            to_string(i) +  "_" + to_string(p) + ".res";
                    nFDeadDens    = "../OutputFiles/deadDens_" + to_string(i) +
                            "_" + to_string(p) + ".res";
                    nFCycle       = "../OutputFiles/cycle_" + to_string(i) +
                            "_" + to_string(p) + ".res";
                    nFHypDens     = "../OutputFiles/hypDens_" + to_string(i) +
                            "_" + to_string(p) + ".res";
                    nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(i) +
                            "_" + to_string(p) + ".res";
                    nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(i) +
                            "_" + to_string(p) + ".res";
                    createInFiles(nrow, ncol, nlayer, tumDens, sigmaTum,
                                  vascDens, inTum, inVes);
                    model(x, y[p], nrow, ncol, nlayer, cellSize, inTum, inVes,
                          nFTumDens, nFTumVol, nFVascDens, nFKilledCells,
                          nFDeadDens, nFCycle, nFHypDens, nFPO2Stat,
                          nFVegfStat);
                    nEv++;

                    cout << nEv << " out of " << nEvTot <<
                            " evaluations of the model" << endl;
                    cout << "-----------------------------------------" << endl;


                    for(int j(0); j < nOut; j++){
                        ymean[j] += y[p][j];
                    }
                }

                for(int j(0); j < nOut; j++){
                    ymean[j] /= P;
                }

                f8wTumDens     << ymean[0]   << " ";
                f12wTumDens    << ymean[1]   << " ";
                f8wTumVol      << ymean[2]   << " ";
                f12wTumVol     << ymean[3]   << " ";
                f8wIntTumDens  << ymean[4]   << " ";
                f12wIntTumDens << ymean[5]   << " ";
                f8wIntTumVol   << ymean[6]   << " ";
                f12wIntTumVol  << ymean[7]   << " ";
                fKilled50      << ymean[8]   << " ";
                fKilled80      << ymean[9]   << " ";
                fKilled90      << ymean[10]  << " ";
                fKilled95      << ymean[11]  << " ";
                fTimeTo95      << ymean[12]  << " ";
                fKilled99      << ymean[13]  << " ";
                fTimeTo99      << ymean[14]  << " ";
                fKilled999     << ymean[15]  << " ";
                fRec           << ymean[16]  << " ";
                fRecTumDens    << ymean[17]  << " ";
                fRecTumVol     << ymean[18]  << " ";
                fRecTime       << ymean[19]  << " ";

                if(P > 1){
                    for(int j(0); j < nOut; j++){
                        for(int p(0); p < P; p++){
                            ystd[j] += (y[p][j] - ymean[j]) * (y[p][j] -
                                                               ymean[j]);
                        }
                        ystd[j] = sqrt(ystd[j] / (P - 1.0));
                    }

                    f8wTumDens     << ystd[0]   << " ";
                    f12wTumDens    << ystd[1]   << " ";
                    f8wTumVol      << ystd[2]   << " ";
                    f12wTumVol     << ystd[3]   << " ";
                    f8wIntTumDens  << ystd[4]   << " ";
                    f12wIntTumDens << ystd[5]   << " ";
                    f8wIntTumVol   << ystd[6]   << " ";
                    f12wIntTumVol  << ystd[7]   << " ";
                    fKilled50      << ystd[8]   << " ";
                    fKilled80      << ystd[9]   << " ";
                    fKilled90      << ystd[10]  << " ";
                    fKilled95      << ystd[11]  << " ";
                    fTimeTo95      << ystd[12]  << " ";
                    fKilled99      << ystd[13]  << " ";
                    fTimeTo99      << ystd[14]  << " ";
                    fKilled999     << ystd[15]  << " ";
                    fRec           << ystd[16]  << " ";
                    fRecTumDens    << ystd[17]  << " ";
                    fRecTumVol     << ystd[18]  << " ";
                    fRecTime       << ystd[19]  << " ";
                }

                f8wTumDens     << endl;
                f12wTumDens    << endl;
                f8wTumVol      << endl;
                f12wTumVol     << endl;
                f8wIntTumDens  << endl;
                f12wIntTumDens << endl;
                f8wIntTumVol   << endl;
                f12wIntTumVol  << endl;
                fKilled50      << endl;
                fKilled80      << endl;
                fKilled90      << endl;
                fKilled95      << endl;
                fTimeTo95      << endl;
                fKilled99      << endl;
                fTimeTo99      << endl;
                fKilled999     << endl;
                fRec           << endl;
                fRecTumDens    << endl;
                fRecTumVol     << endl;
                fRecTime       << endl;

                fCombDens << tumDens << " " <<  sigmaTum << " " << vascDens <<
                             endl;
                vascDens += hVascDens;
                i++;
            }
            sigmaTum += hSigmaTum;
        }
        tumDens += hTumDens;
    }

    f8wTumDens.close();
    f12wTumDens.close();
    f8wTumVol.close();
    f12wTumVol.close();
    f8wIntTumDens.close();
    f12wIntTumDens.close();
    f8wIntTumVol.close();
    f12wIntTumVol.close();
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
    fRecTumVol.close();
    fRecTime.close(););
    fRecTime.close();
    fCombDens.close();
}


/*------------------------------------------------------------------------------
 * This function studies the impact of varying one or more parameters on both
 * scalar and time-dependent outputs of the model. For the most relevant
 * parameters of the model, several values within their ranges are used. For the
 * least relevant parameters, fixed values are considered.
 *
 * Inputs:
 *  - nFVarPar: name of the file containing a list of the varying parameters and
 *  their possible values,
 *  - nFMostRelPar: name of the file containing a list of the the most relevant
 *  parameters and their ranges,
 *  - nFLeastRelPar: name of the file containing a list of the the least
 *  relevant parameters and their values,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration,
 *  - L: number of possible values of the most relevant parameters,
 *  - P: number of repetitions for each combination of parameters.
------------------------------------------------------------------------------*/

void varErr(const string nFVarPar, const string nFMostRelPar,
            const string nFLeastPar, const string nFInTissueDim,
            const string nFInTum, const string nFInVes,
            const int L, const int P){
    const int K(39), nOut(20);
    int nLeastRelPar, nMostRelPar, nVarPar;

    ifstream fMostRelPar(nFMostRelPar.c_str());

    fMostRelPar >> nMostRelPar;

    int Xi[nMostRelPar];
    double h, x[K], X[nMostRelPar][L];

    for(int k(0); k < nMostRelPar; k++){
        fMostRelPar >> Xi[k] >> X[k][0] >> X[k][L - 1];
        h = (X[k][L - 1] - X[k][0]) / (L - 1);
        for(int l(1); l < L - 1; l++){
            X[k][l] = X[k][0] + l * h;
        }
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

    const int nEvTot(2 * pow(L, nMostRelPar) * P);
    int nEv(0);
    double y0[P][nOut], y0mean[nOut], y0std[nOut];
    double y1[P][nOut], y1mean[nOut], y1std[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream f8wTumDens("../OutputFiles/8wTumDens.res");
    ofstream f12wTumDens("../OutputFiles/12wTumDens.res");
    ofstream f8wTumVol("../OutputFiles/8wTumVol.res");
    ofstream f12wTumVol("../OutputFiles/12wTumVol.res");
    ofstream f8wIntTumDens("../OutputFiles/8wIntTumDens.res");
    ofstream f12wIntTumDens("../OutputFiles/12wIntTumDens.res");
    ofstream f8wIntTumVol("../OutputFiles/8wIntTumVol.res");
    ofstream f12wIntTumVol("../OutputFiles/12wIntTumVol.res");
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
    ofstream fRecTumVol("../OutputFiles/recTumVol.res");
    ofstream fRecTime("../OutputFiles/recTime.res");
    ofstream fCombPar("../OutputFiles/combPar.res");

    int count(0);

    for(int i1(0); i1 < L; i1++){
        x[Xi[0]] = X[0][i1];
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
            nFVascDens    = "../OutputFiles/vascDens_" + to_string(count) +
                    "_" + "0_" + to_string(p) + ".res";
            nFKilledCells = "../OutputFiles/killedCells_" + to_string(count) +
                    "_" + "0_" + to_string(p) + ".res";
            nFDeadDens    = "../OutputFiles/deadDens_" + to_string(count) +
                    "_" + "0_" + to_string(p) + ".res";
            nFCycle       = "../OutputFiles/cycle_" + to_string(count) + "_" +
                    "0_" + to_string(p) + ".res";
            nFHypDens     = "../OutputFiles/hypDens_" + to_string(count) + "_" +
                    "0_" + to_string(p) + ".res";
            nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(count) + "_" +
                    "0_" + to_string(p) + ".res";
            nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(count) +
                    "_" + "0_" + to_string(p) + ".res";

            model(x, y0[p], nrow, ncol, nlayer, cellSize, inTum, inVes,
                  nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens,
                  nFCycle, nFHypDens, nFPO2Stat, nFVegfStat);
            nEv++;

            cout << nEv << " out of " << nEvTot <<
                    " evaluations of the model" << endl;
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
            nFVascDens    = "../OutputFiles/vascDens_" + to_string(count) +
                    "_" + "1_" + to_string(p) + ".res";
            nFKilledCells = "../OutputFiles/killedCells_" + to_string(count) +
                    "_" + "1_" + to_string(p) + ".res";
            nFDeadDens    = "../OutputFiles/deadDens_" + to_string(count) +
                    "_" + "1_" + to_string(p) + ".res";
            nFCycle       = "../OutputFiles/cycle_" + to_string(count) + "_" +
                    "1_" + to_string(p) + ".res";
            nFHypDens     = "../OutputFiles/hypDens_" + to_string(count) + "_" +
                    "1_" + to_string(p) + ".res";
            nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(count) + "_" +
                    "1_" + to_string(p) + ".res";
            nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(count) +
                    "_" + "1_" + to_string(p) + ".res";

            model(x, y1[p], nrow, ncol, nlayer, cellSize, inTum, inVes,
                  nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens,
                  nFCycle, nFHypDens, nFPO2Stat, nFVegfStat);
            nEv++;

            cout << nEv << " out of " << nEvTot <<
                    " evaluations of the model" << endl;
            cout << "---------------------------------------------" << endl;

            for(int j(0); j < nOut; j++){
                y1mean[j] += y1[p][j];
            }
        }
        count++;

        for(int j(0); j < nOut; j++){
            y1mean[j] /= P;
        }

        f8wTumDens     << y0mean[0]   << " " << y1mean[0]    << " ";
        f12wTumDens    << y0mean[1]   << " " << y1mean[1]    << " ";
        f8wTumVol      << y0mean[2]   << " " << y1mean[2]    << " ";
        f12wTumVol     << y0mean[3]   << " " << y1mean[3]    << " ";
        f8wIntTumDens  << y0mean[4]   << " " << y1mean[4]    << " ";
        f12wIntTumDens << y0mean[5]   << " " << y1mean[5]    << " ";
        f8wIntTumVol   << y0mean[6]   << " " << y1mean[6]    << " ";
        f12wIntTumVol  << y0mean[7]   << " " << y1mean[7]    << " ";
        fKilled50      << y0mean[8]   << " " << y1mean[8]    << " ";
        fKilled80      << y0mean[9]   << " " << y1mean[9]    << " ";
        fKilled90      << y0mean[10]  << " " << y1mean[10]   << " ";
        fKilled95      << y0mean[11]  << " " << y1mean[11]   << " ";
        fTimeTo95      << y0mean[12]  << " " << y1mean[12]   << " ";
        fKilled99      << y0mean[13]  << " " << y1mean[13]   << " ";
        fTimeTo99      << y0mean[14]  << " " << y1mean[14]   << " ";
        fKilled999     << y0mean[15]  << " " << y1mean[15]   << " ";
        fRec           << y0mean[16]  << " " << y1mean[16]   << " ";
        fRecTumDens    << y0mean[17]  << " " << y1mean[17]   << " ";
        fRecTumVol     << y0mean[18]  << " " << y1mean[18]   << " ";
        fRecTime       << y0mean[19]  << " " << y1mean[19]   << " ";


        if(P > 1){
            for(int j(0); j < nOut; j++){
                for(int p(0); p < P; p++){
                    y0std[j] += (y0[p][j] - y0mean[j]) * (y0[p][j] - y0mean[j]);
                    y1std[j] += (y1[p][j] - y1mean[j]) * (y1[p][j] - y1mean[j]);
                }
                y0std[j] = sqrt(y0std[j] / (P - 1.0));
                y1std[j] = sqrt(y1std[j] / (P - 1.0));
            }

            f8wTumDens     << y0std[0]   << " " << y1std[0]    << " ";
            f12wTumDens    << y0std[1]   << " " << y1std[1]    << " ";
            f8wTumVol      << y0std[2]   << " " << y1std[2]    << " ";
            f12wTumVol     << y0std[3]   << " " << y1std[3]    << " ";
            f8wIntTumDens  << y0std[4]   << " " << y1std[4]    << " ";
            f12wIntTumDens << y0std[5]   << " " << y1std[5]    << " ";
            f8wIntTumVol   << y0std[6]   << " " << y1std[6]    << " ";
            f12wIntTumVol  << y0std[7]   << " " << y1std[7]    << " ";
            fKilled50      << y0std[8]   << " " << y1std[8]    << " ";
            fKilled80      << y0std[9]   << " " << y1std[9]    << " ";
            fKilled90      << y0std[10]  << " " << y1std[10]   << " ";
            fKilled95      << y0std[11]  << " " << y1std[11]   << " ";
            fTimeTo95      << y0std[12]  << " " << y1std[12]   << " ";
            fKilled99      << y0std[13]  << " " << y1std[13]   << " ";
            fTimeTo99      << y0std[14]  << " " << y1std[14]   << " ";
            fKilled999     << y0std[15]  << " " << y1std[15]   << " ";
            fRec           << y0std[16]  << " " << y1std[16]   << " ";
            fRecTumDens    << y0std[17]  << " " << y1std[17]   << " ";
            fRecTumVol     << y0std[18]  << " " << y1std[18]   << " ";
            fRecTime       << y0std[19]  << " " << y1std[19]   << " ";
        }

        f8wTumDens     << endl;
        f12wTumDens    << endl;
        f8wTumVol      << endl;
        f12wTumVol     << endl;
        f8wIntTumDens  << endl;
        f12wIntTumDens << endl;
        f8wIntTumVol   << endl;
        f12wIntTumVol  << endl;
        fKilled50      << endl;
        fKilled80      << endl;
        fKilled90      << endl;
        fKilled95      << endl;
        fTimeTo95      << endl;
        fKilled99      << endl;
        fTimeTo99      << endl;
        fKilled999     << endl;
        fRec           << endl;
        fRecTumDens    << endl;
        fRecTumVol     << endl;
        fRecTime       << endl;

        for (int j(0); j < nMostRelPar; j ++){
            fCombPar << x[Xi[j]] << " ";
        }
        fCombPar << endl;
    }

    f8wTumDens.close();
    f12wTumDens.close();
    f8wTumVol.close();
    f12wTumVol.close();
    f8wIntTumDens.close();
    f12wIntTumDens.close();
    f8wIntTumVol.close();
    f12wIntTumVol.close();
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
    fRecTumVol.close();
    fRecTime.close(););
    fRecTime.close();
    fCombPar.close();
}


/*------------------------------------------------------------------------------
 * This function compares the scalar output values obtained for two or more
 * evaluations of the model using parameters defined in input files.
 *
 * Inputs:
 *  - nFPar: vector with the names of the files containing the values of the
 *  parameters of the model,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration.
------------------------------------------------------------------------------*/

void varParFromFiles(const vector<string> nFPar, const string nFInTissueDim,
                     const string nFInTum, const string nFInVes){
    const int K(38), L(nFPar.size()), nOut(20);
    double x[K], y[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;
    ofstream f8wTumDens, f12wTumDens, f8wTumVol, f12wTumVol;
    ofstream f8wIntTumDens, f12wIntTumDens, f8wIntTumVol, f12wIntTumVol;
    ofstream fKilled50, fKilled80, fKilled90, fKilled95, fKilled99, fKilled999;
    ofstream fTimeTo95, fTimeTo99;
    ofstream fRec, fRecTumDens, fRecTumVol, fRecTime;

    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer, cellSize,
                inTum, inVes);

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
        nFDeadDens    = "../OutputFiles/deadDens_" + to_string(i) + ".res";
        nFCycle       = "../OutputFiles/cycle_" + to_string(i) + ".res";
        nFHypDens     = "../OutputFiles/hypDens_" + to_string(i) + ".res";
        nFPO2Stat     = "../OutputFiles/pO2Stat_" + to_string(i) + ".res";
        nFVegfStat    = "../OutputFiles/vegfStat_" + to_string(i) + ".res";

        f8wTumDens.open("../OutputFiles/8wTumDens_" + to_string(i) +
                              ".res");
        f12wTumDens.open("../OutputFiles/12wTumDens_" + to_string(i) +
                          ".res");
        f8wTumVol.open("../OutputFiles/8wTumVol_" + to_string(i) +
                              ".res");
        f12wTumVol.open("../OutputFiles/12wTumVol_" + to_string(i) +
                          ".res");
        f8wIntTumDens.open("../OutputFiles/8wIntTumDens_" + to_string(i) +
                              ".res");
        f12wIntTumDens.open("../OutputFiles/12wIntTumDens_" + to_string(i) +
                          ".res");
        f8wIntTumVol.open("../OutputFiles/8wIntTumVol_" + to_string(i) +
                              ".res");
        f12wIntTumVol.open("../OutputFiles/12wIntTumVol_" + to_string(i) +
                          ".res");
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
        fRecTumVol.open("../OutputFiles/recTumVol_" + to_string(i) + ".res");
        fRecTime.open("../OutputFiles/recTime_" + to_string(i) + ".res");

        model(x, y, nrow, ncol, nlayer, cellSize, inTum, inVes, nFTumDens,
              nFTumVol, nFVascDens, nFKilledCells, nFDeadDens, nFCycle,
              nFHypDens, nFPO2Stat, nFVegfStat);

        f8wTumDens     << y[0];
        f12wTumDens    << y[1];
        f8wTumVol      << y[2];
        f12wTumVol     << y[3];
        f8wIntTumDens  << y[4];
        f12wIntTumDens << y[5];
        f8wIntTumVol   << y[6];
        f12wIntTumVol  << y[7];
        fKilled50      << y[8];
        fKilled80      << y[9];
        fKilled90      << y[10];
        fKilled95      << y[11];
        fTimeTo95      << y[12];
        fKilled99      << y[13];
        fTimeTo99      << y[14];
        fKilled999     << y[15];
        fRec           << y[16];
        fRecTumDens    << y[17];
        fRecTumVol     << y[18];
        fRecTime       << y[19];

        f8wTumDens.close();
        f12wTumDens.close();
        f8wTumVol.close();
        f12wTumVol.close();
        f8wIntTumDens.close();
        f12wIntTumDens.close();
        f8wIntTumVol.close();
        f12wIntTumVol.close();
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
        fRecTumVol.close();
        fRecTime.close(););
        fRecTime.close();

        cout << i + 1 << " out of " << L << " evaluations of the model" << endl;
        cout << "---------------------------------------------" << endl;
    }
}


/*------------------------------------------------------------------------------
 * This function compares the scalar output values obtained for two or more
 * evaluations of the model using parameters defined in input files.
 *
 * Inputs:
 *  - nFPar: vector with the names of the files containing the values of the
 *  parameters of the model,
 *  - nFInTissuePar: name of the file containing the parameters of an
 *  artificial tissue,
------------------------------------------------------------------------------*/

void varParFromFiles(const string nFInTissuePar, const string nFPar,
                     const string nFTreatment){
    const int K(38), nOut(20);
    double x[K], y[nOut];
    string nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens;
    string nFCycle, nFHypDens, nFPO2Stat, nFVegfStat;

    int nrow, ncol, nlayer;
    double cellSize, tumArea, tumDens, vascDens;
    vector<bool> inTum, inVes;
    Treatment treatment;

    readInFiles(nFInTissuePar, nFTreatment, cellSize, tumArea, tumDens,
                vascDens, treatment);
    createInFiles(cellSize, tumArea, tumDens, vascDens, nrow, ncol, nlayer,
                  inTum, inVes);

    cout << "cell size: " << cellSize << endl;
    cout << "tumor area: " << tumArea << endl;
    cout << "tumor density in the tumor area: " << tumDens << endl;
    cout << "vascular density: " << vascDens << endl;
    cout << "nrow: " << nrow << endl;
    cout << "ncol: " << ncol << endl;
    cout << "nlayer: " << nlayer << endl;

    ifstream fPar(nFPar.c_str());
    for(int k(0); k < K; k++){
        fPar >> x[k];
    }
    fPar.close();

    nFTumDens     = "../OutputFiles/tumDens.res";
    nFTumVol      = "../OutputFiles/tumVol.res";
    nFVascDens    = "../OutputFiles/vascDens.res";
    nFKilledCells = "../OutputFiles/killedCells.res";
    nFDeadDens    = "../OutputFiles/deadDens.res";
    nFCycle       = "../OutputFiles/cycle.res";
    nFHypDens     = "../OutputFiles/hypDens.res";
    nFPO2Stat     = "../OutputFiles/pO2Stat.res";
    nFVegfStat    = "../OutputFiles/vegfStat.res";

    model(x, y, nrow, ncol, nlayer, cellSize, inTum, inVes, &treatment,
          nFTumDens, nFTumVol, nFVascDens, nFKilledCells, nFDeadDens, nFCycle,
          nFHypDens, nFPO2Stat, nFVegfStat);



    ofstream f8wTumDens("../OutputFiles/8wTumDens.res");
    ofstream f12wTumDens("../OutputFiles/12wTumDens.res");
    ofstream f8wTumVol("../OutputFiles/8wTumVol.res");
    ofstream f12wTumVol("../OutputFiles/12wTumVol.res");
    ofstream f8wIntTumDens("../OutputFiles/8wIntTumDens.res");
    ofstream f12wIntTumDens("../OutputFiles/12wIntTumDens.res");
    ofstream f8wIntTumVol("../OutputFiles/8wIntTumVol.res");
    ofstream f12wIntTumVol("../OutputFiles/12wIntTumVol.res");
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
    ofstream fRecTumVol("../OutputFiles/recTumVol.res");
    ofstream fRecTime("../OutputFiles/recTime.res");

    f8wTumDens     << y[0];
    f12wTumDens    << y[1];
    f8wTumVol      << y[2];
    f12wTumVol     << y[3];
    f8wIntTumDens  << y[4];
    f12wIntTumDens << y[5];
    f8wIntTumVol   << y[6];
    f12wIntTumVol  << y[7];
    fKilled50      << y[8];
    fKilled80      << y[9];
    fKilled90      << y[10];
    fKilled95      << y[11];
    fTimeTo95      << y[12];
    fKilled99      << y[13];
    fTimeTo99      << y[14];
    fKilled999     << y[15];
    fRec           << y[16];
    fRecTumDens    << y[17];
    fRecTumVol     << y[18];
    fRecTime       << y[19];

    f8wTumDens.close();
    f12wTumDens.close();
    f8wTumVol.close();
    f12wTumVol.close();
    f8wIntTumDens.close();
    f12wIntTumDens.close();
    f8wIntTumVol.close();
    f12wIntTumVol.close();
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
    fRecTumVol.close();
    fRecTime.close();
    cout << "---------------------------------------------" << endl;
}

/*------------------------------------------------------------------------------
 * This function evaluates the model for random combinations of the values of
 * the parameters within their ranges. The scalar ouptuts are calculated.
 *
 * Inputs:
 *  - N: number of random combinations,
 *  - P: number of repetitions for each combination,
 *  - nRefParInt: name of the file containing the reference ranges for all the
 *  parameters,
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration.
------------------------------------------------------------------------------*/

void varStoch(const int N, const int P, const string nFRefParInt,
              const string nFInTissueDim, const string nFInTum,
              const string nFInVes){
    int nrow, ncol, nlayer;
    double cellSize;
    vector<bool> inTum, inVes;

    readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                cellSize, inTum, inVes);

    const int K(39), nOut(20);
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

    ofstream f8wTumDens("../OutputFiles/8wTumDens.res");
    ofstream f12wTumDens("../OutputFiles/12wTumDens.res");
    ofstream f8wTumVol("../OutputFiles/8wTumVol.res");
    ofstream f12wTumVol("../OutputFiles/12wTumVol.res");
    ofstream f8wIntTumDens("../OutputFiles/8wIntTumDens.res");
    ofstream f12wIntTumDens("../OutputFiles/12wIntTumDens.res");
    ofstream f8wIntTumVol("../OutputFiles/8wIntTumVol.res");
    ofstream f12wIntTumVol("../OutputFiles/12wIntTumVol.res");
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
    ofstream fRecTumVol("../OutputFiles/recTumVol.res");
    ofstream fRecTime("../OutputFiles/recTime.res");

    for(int i(0); i < N; i++){
        for(int l(0); l < P; l++){
            model(X[i], Y[l][i], nrow, ncol, nlayer, cellSize, inTum, inVes);
            nEv++;
            cout << nEv << " out of " << nEvTot <<
                    " evaluations of the model" << endl;
            cout << "---------------------------------------------" << endl;

            f8wTumDens     << Y[l][i][0]  << " ";
            f12wTumDens    << Y[l][i][1]  << " ";
            f8wTumVol      << Y[l][i][2]  << " ";
            f12wTumVol     << Y[l][i][3]  << " ";
            f8wIntTumDens  << Y[l][i][4]  << " ";
            f12wIntTumDens << Y[l][i][5]  << " ";
            f8wIntTumVol   << Y[l][i][6]  << " ";
            f12wIntTumVol  << Y[l][i][7]  << " ";
            fKilled50      << Y[l][i][8]  << " ";
            fKilled80      << Y[l][i][9]  << " ";
            fKilled90      << Y[l][i][10] << " ";
            fKilled95      << Y[l][i][11] << " ";
            fTimeTo95      << Y[l][i][12] << " ";
            fKilled99      << Y[l][i][13] << " ";
            fTimeTo99      << Y[l][i][14] << " ";
            fKilled999     << Y[l][i][15] << " ";
            fRec           << Y[l][i][16] << " ";
            fRecTumDens    << Y[l][i][17] << " ";
            fRecTumVol     << Y[l][i][18] << " ";
            fRecTime       << Y[l][i][19] << " ";
        }

        f8wTumDens     << endl;
        f12wTumDens    << endl;
        f8wTumVol      << endl;
        f12wTumVol     << endl;
        f8wIntTumDens  << endl;
        f12wIntTumDens << endl;
        f8wIntTumVol   << endl;
        f12wIntTumVol  << endl;
        fKilled50      << endl;
        fKilled80      << endl;
        fKilled90      << endl;
        fKilled95      << endl;
        fTimeTo95      << endl;
        fKilled99      << endl;
        fTimeTo99      << endl;
        fKilled999     << endl;
        fRec           << endl;
        fRecTumDens    << endl;
        fRecTumVol     << endl;
        fRecTime       << endl;
    }

    f8wTumDens.close();
    f12wTumDens.close();
    f8wTumVol.close();
    f12wTumVol.close();
    f8wIntTumDens.close();
    f12wIntTumDens.close();
    f8wIntTumVol.close();
    f12wIntTumVol.close();
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
    fRecTumVol.close();
    fRecTime.close(););
    fRecTime.close();

    free2D(X, N);
    free3D(Y, P, N);
}

