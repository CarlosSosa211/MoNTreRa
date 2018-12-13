#include "sensAn.hpp"

using namespace std;

void evalR(const int nMethod, const int nModel){
    string nFDim;

    switch(nMethod){
    case 0:
        nFDim = "../InputFiles/morrisDim.dat";
        break;
    case 1:
        nFDim = "../InputFiles/sobolDim.dat";
    }

    int K;
    double doubN;
    ifstream fDim(nFDim.c_str());

    fDim >> K >> doubN;
    fDim.close();

    int N(doubN), nEv;

    switch(nMethod){
    case 0:
        nEv = (K + 1) * N;
        break;
    case 1:
        nEv = (K + 2) * N;
        break;
    }

    double x[K];
    ifstream fX("../InputFiles/X.dat");

    switch(nModel){
    case 0:{
        const int nOut(1);
        double y[nOut];
        ofstream fY("../OutputFiles/Y.res");

        for(int i(0); i < nEv; i++){
            for(int k(0); k < K; k++){
                fX >> x[k];
            }
            toyModel(x, y);
            //cout << i + 1 << " out of " << nEv << " evaluations of the model" << endl;
            //cout << "---------------------------------------------" << endl;

            fY << y[0] << endl;
        }

        fY.close();
        break;
    }

    case 1:{
        const int nOut(8);
        double y[nOut];
        ofstream fYEndTreatTumDens("../OutputFiles/YEndTreatTumDens.res");
        ofstream fY3MonTumDens("../OutputFiles/Y3MonTumDens.res");
        ofstream fYRecTumDens("../OutputFiles/YRecTumDens.res");
        ofstream fYTumVol("../OutputFiles/YTumVol.res");
        ofstream fYIntTumDens("../OutputFiles/YIntTumDens.res");
        ofstream fYTimeTo95("../OutputFiles/YTimeTo95.res");
        ofstream fYTimeTo99("../OutputFiles/YTimeTo99.res");
        ofstream fYRecTime("../OutputFiles/YRecTime.res");

        for(int i(0); i < nEv; i++){
            for(int k(0); k < K; k++){
                fX >> x[k];
            }
            model(x, y);
            //cout << i + 1 << " out of " << nEv << " evaluations of the model" << endl;
            //cout << "---------------------------------------------" << endl;

            fYEndTreatTumDens << y[0] << endl;
            fY3MonTumDens     << y[1] << endl;
            fYRecTumDens      << y[2] << endl;
            fYTumVol          << y[3] << endl;
            fYIntTumDens      << y[4] << endl;
            fYTimeTo95        << y[5] << endl;
            fYTimeTo99        << y[6] << endl;
            fYRecTime         << y[7] << endl;
        }

        fYEndTreatTumDens.close();
        fY3MonTumDens.close();
        fYRecTumDens.close();
        fYTumVol.close();
        fYIntTumDens.close();
        fYTimeTo95.close();
        fYTimeTo99.close();
        fYRecTime.close();

        break;
    }
    }

    fX.close();
}
