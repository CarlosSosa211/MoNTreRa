#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>

#include "morris.hpp"
#include "sensAnR.hpp"
#include "sobol.hpp"
#include "tcp.hpp"
#include "var.hpp"

using namespace std;

int main(){
    //const int N(100);
    const int kp(5), L(5), p(20), P(5), N(10);
    //const int nMethod(0), nModel(1);
    //string nFRefParInt("../InputFiles/refParIntOneAlphaBeta.dat");
    //string nFRefParInt("../InputFiles/refParIntRT.dat");
    //string nFRefParInt("../InputFiles/refParIntToy.dat");
    //string nRefParMean("../InputFiles/refParMeanRT.dat");
    string nFMostRelPar("../InputFiles/mostRelParArrest.dat");
    string nFLeastRelPar("../InputFiles/leastRelParArrest.dat");
    string nFVarPar("../InputFiles/varParArrest.dat");
    string nFInTissueDim("../InputFiles/tissueDim.dat");
    string nFInTum("../InputFiles/inTum.dat");
    string nFInVes("../InputFiles/inVes.dat");
    string nFInPO2("../InputFiles/inPO2.dat");
    //vector<string> nFPar;
    //nFPar.push_back("../InputFiles/par37_2.dat");
    //nFPar.push_back("../InputFiles/par20_3.dat");
    /*string nFInTissueTCP("../InputFiles/inTissueTCP0.dat");
    string nFParTCP("../InputFiles/parTCP.dat");
    vector<string> nFTreatmentTCP;
    nFTreatmentTCP.push_back("../InputFiles/1MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/2MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/3MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/4MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/5MonFri.dat");*/

    srand(time(NULL));
    //evalR(nMethod, nModel, nFInTissueDim, nFInTum, nFInVes, nFInPO2);
    //morrisRT(N, p, nFRefParInt, nFInTissueDim, nFInTum, nFInVes, nFInPO2);
    //morrisToy(N, p, nFRefParInt);
    /*morrisVarRangeRT(kp, L, N, p, nFRefParInt, nFInTissueDim, nFInTum,
        nFInVes, nFInPO2);*/
    //morrisVarRangeToy(kp, L, N, p, nFRefParInt);
    //sobolRT(N, nFRefParInt, nFInTissueDim, nFInTum, nFInVes, nFInPO2);
    //sobolToy(N, nFRefParInt);
    //sobolFromFiles(2);
    //tcp(N, nFInTissueTCP, nFParTCP, nFTreatmentTCP);
    /*tcp(N, nFInTissueTCP, nFParTCP, nFTreatmentTCP, nFInTissueDim, nFInTum,
        nFInVes, nFInPO2);*/
    /*var1ParRange(kp, L, nFRefParInt, nFInTissueDim, nFInTum, nFInVes,
    nFInPO2);*/
    varErr(nFVarPar, nFMostRelPar, nFLeastRelPar, nFInTissueDim, nFInTum,
           nFInVes, nFInPO2, L, P);
    //varParFromFiles(nFPar, nFInTissueDim, nFInTum, nFInVes, nFInPO2);
    //varStoch(N, P, nFRefParInt, nFInTissueDim, nFInTum, nFInVes, nFInPO2);
}



