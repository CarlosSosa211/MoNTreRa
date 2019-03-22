#include "sensAn.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * This functions evaluates a model using the combinations of parameters already
 * built in R and writes the outputs in files.
 *
 * Inputs:
 *  - nMethod: method used (0 for Morris and 1 for Sobol),
 *  - nModel: model considered (0 for the model of tumour growth and response to
 *  radiotherapy and 1 for the toy model),
 *  - nFInTissueDim: name of the file containing the dimensions of a
 *  histological specimen; it is empty if the toy model is considered,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration; it is empty if the toy model is considered,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration; it is empty if the toy model is considered.
------------------------------------------------------------------------------*/

void evalR(const int nMethod, const int nModel, const std::string nFInTissueDim,
           const std::string nFInTum, const std::string nFInVes){
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
            cout << i + 1 << " out of " << nEv << " evaluations of the model" <<
                    endl;
            cout << "---------------------------------------------" << endl;

            fY << y[0] << endl;
        }

        fY.close();
        break;
    }

    case 1:{

        int nrow, ncol, nlayer;
        double cellSize;
        vector<bool> inTum, inVes;

        if(!nFInTissueDim.empty() && !nFInTum.empty() && !nFInVes.empty()){
            readInFiles(nFInTissueDim, nFInTum, nFInVes, nrow, ncol, nlayer,
                        cellSize, inTum, inVes);
        }

        const int nOut(15);
        double y[nOut];
        ofstream fYEndTreatTumDens("../OutputFiles/YEndTreatTumDens.res");
        ofstream fY3MonTumDens("../OutputFiles/Y3MonTumDens.res");
        ofstream fYTumVol("../OutputFiles/YTumVol.res");
        ofstream fYIntTumDens("../OutputFiles/YIntTumDens.res");
        ofstream fYKilled50("../OutputFiles/YKilled50.res");
        ofstream fYKilled80("../OutputFiles/YKilled80.res");
        ofstream fYKilled90("../OutputFiles/YKilled90.res");
        ofstream fYKilled95("../OutputFiles/YKilled95.res");
        ofstream fYTimeTo95("../OutputFiles/YTimeTo95.res");
        ofstream fYKilled99("../OutputFiles/YKilled99.res");
        ofstream fYTimeTo99("../OutputFiles/YTimeTo99.res");
        ofstream fYKilled999("../OutputFiles/YKilled999.res");
        ofstream fYRec("../OutputFiles/YRec.res");
        ofstream fYRecTumDens("../OutputFiles/YRecTumDens.res");
        ofstream fYRecTime("../OutputFiles/YRecTime.res");

        for(int i(0); i < nEv; i++){
            for(int k(0); k < K; k++){
                fX >> x[k];
            }
            model(x, y, nrow, ncol, nlayer, cellSize, inTum, inVes);
            cout << i + 1 << " out of " << nEv << " evaluations of the model" <<
                    endl;
            cout << "---------------------------------------------" << endl;

            fYEndTreatTumDens << y[0]  << endl;
            fY3MonTumDens     << y[1]  << endl;
            fYTumVol          << y[2]  << endl;
            fYIntTumDens      << y[3]  << endl;
            fYKilled50        << y[4]  << endl;
            fYKilled80        << y[5]  << endl;
            fYKilled90        << y[6]  << endl;
            fYKilled95        << y[7]  << endl;
            fYTimeTo95        << y[8]  << endl;
            fYKilled99        << y[9]  << endl;
            fYTimeTo99        << y[10] << endl;
            fYKilled999       << y[11] << endl;
            fYRec             << y[12] << endl;
            fYRecTumDens      << y[13] << endl;
            fYRecTime         << y[14] << endl;
        }

        fYEndTreatTumDens.close();
        fY3MonTumDens.close();
        fYTumVol.close();
        fYIntTumDens.close();
        fYKilled50.close();
        fYKilled80.close();
        fYKilled90.close();
        fYKilled95.close();
        fYTimeTo95.close();
        fYKilled99.close();
        fYTimeTo99.close();
        fYKilled999.close();
        fYRec.close();
        fYRecTumDens.close();
        fYRecTime.close();

        break;
    }
    }

    fX.close();
}
