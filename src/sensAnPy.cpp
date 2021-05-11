#include "evalModel.hpp"
#include "sensAnPy.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * This functions evaluates a model using the combinations of parameters already
 * built in Python and writes the outputs in files.
 *
 * Inputs:
 *  - nMethod: method used (0 for Morris and 1 for Sobol),
 *  - nFX: name of the file containing the model parameter values,
 *  - nFY: name of the file containing the model output valus.
------------------------------------------------------------------------------*/

void evalPy(const int nMethod, const string nFX, const string nFY){
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
    ifstream fX(nFX);

    const int nOut(20);
    double y[nOut];
    ofstream fY(nFY.c_str());

    for(int i(0); i < nEv; i++){
        for(int k(0); k < K; k++){
            fX >> x[k];
        }
        reducedArtModel(x, y);
        cout << i + 1 << " out of " << nEv << " evaluations of the model" <<
                endl;
        cout << "---------------------------------------------" << endl;

        for(int j(0); j < nOut; j++){
            fY << y[j] << " ";
        }
        fY << endl;
    }

    fX.close();
    fY.close();
}

/*------------------------------------------------------------------------------
 * This functions evaluates a model using the combinations of parameters already
 * built in Python and writes the outputs in files.
 *
 * Inputs:
 *  - nMethod: method used (0 for Morris and 1 for Sobol),
 *  - nFX: name of the file containing the model parameter values,
 *  - nFInTissuePar: name of the file containing the parameters of an
 *  artificial tissue,
 *  -  nFTreatment: name of the file containing the administered treatment,
 *  - nFY: name of the file containing the model output valus.
------------------------------------------------------------------------------*/

void evalPy(const int nMethod, const string nFInTissuePar, const string nFX,
            const string nFTreatment, const string nFY){
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
    ifstream fX(nFX);

    int nrow, ncol, nlayer;
    double cellSize, tumArea, radRatioTum, tumDens, vascDens;
    vector<bool> inTum, inVes;
    Treatment treatment;

    readInFiles(nFInTissuePar, cellSize, tumArea, radRatioTum,
                tumDens, vascDens, nFTreatment, treatment);
    createInFiles(cellSize, tumArea, radRatioTum, tumDens, vascDens, nrow, ncol,
                  nlayer, inTum, inVes);

    cout << nrow << " " << ncol << " " << nlayer << endl;

    const int nOut(20);
    double y[nOut];
    ofstream fY(nFY.c_str());

    for(int i(0); i < nEv; i++){
        for(int k(0); k < K; k++){
            fX >> x[k];
        }
        reducedModel(x, y, nrow, ncol, nlayer, cellSize, inTum, inVes,
                     &treatment);
        cout << i + 1 << " out of " << nEv << " evaluations of the model" <<
                endl;
        cout << "---------------------------------------------" << endl;

        for(int j(0); j < nOut; j++){
            fY << y[j] << " ";
        }
        fY << endl;
    }

    fX.close();
    fY.close();
}
