#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

#include "morris.hpp"
#include "sensAnR.hpp"
#include "sobol.hpp"
#include "tcp.hpp"
#include "var.hpp"

using namespace std;

int main(){
    const int kp(5), L(5), p(20), P(5), N(1000);
    //const int nMethod(0), nModel(1);
    //string nFRefParInt("../InputFiles/refParIntOneAlphaBeta.dat");
    //string nFRefParInt("../InputFiles/refParIntRT.dat");
    //string nFRefParInt("../InputFiles/refParIntRTArt.dat");
    //string nFRefParInt("../InputFiles/refParIntToy.dat");
    //string nFRefParInt("../InputFiles/refParIntRedRT.dat");
    //string nFRefParInt("../InputFiles/refParIntFracRT.dat");
    //string nFRefParMean("../InputFiles/refParMeanRTNoDose.dat");
    /*string nFMostRelPar("../InputFiles/mostRelParAlphaBeta.dat");
    string nFLeastRelPar("../InputFiles/leastRelParAlphaBeta.dat");
    string nFVarPar("../InputFiles/varParAlphaBeta.dat");*/
    //string nFInTissueDim("../InputFiles/tissueDim.dat");
    //string nFArt("../InputFiles/art.dat");
    //string nFDensInt("../InputFiles/densInt.dat");
    //string nFSigmaTumInt("../InputFiles/sigmaTumInt.dat");
    string nFInTissuePar("../InputFiles/tissuePar.dat");
    //string nFInTum("../InputFiles/inTum.dat");
    //string nFInVes("../InputFiles/inVes.dat");
    /*vector<string> nFPar;
    nFPar.push_back(nFRefParMean);*/
    //nFPar.push_back("../InputFiles/par37_2.dat");
    //nFPar.push_back("../InputFiles/par20_3.dat");
    //string nFInTissueTCP("../InputFiles/inTissueTCP0.dat");
    //string nFParTCP("../InputFiles/parTCP.dat");
    //vector<string> nFTreatmentTCP;
    //nFTreatmentTCP.push_back("../InputFiles/1MonFri.dat");
    //nFTreatmentTCP.push_back("../InputFiles/2MonFri.dat");
    /*nFTreatmentTCP.push_back("../InputFiles/3MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/4MonFri.dat");
    nFTreatmentTCP.push_back("../InputFiles/5MonFri.dat");*/
    //string nFPar(nFRefParMean);
    //string nFTreatment("../InputFiles/treatment.dat");

    srand(time(NULL));
    //evalR(nMethod, nModel, nFInTissueDim, nFInTum, nFInVes);
    //morrisRT(N, p, nFRefParInt, nFInTissueDim, nFInTum, nFInVes);
    //morrisFromFiles(N, p);
    //morrisToy(N, p, nFRefParInt);
    /*morrisVarRangeRT(kp, L, N, p, nFRefParInt, nFInTissueDim, nFInTum,
        nFInVes);*/
    //morrisVarRangeToy(kp, L, N, p, nFRefParInt);
    //sobolRT(N, nFRefParInt, nFInTissueDim, nFInTum, nFInVes);
    //sobolToy(N, nFRefParInt);
    //sobolFromFiles();
    //tcp(N, nFInTissueTCP, nFParTCP, nFTreatmentTCP);
    //tcp(N, nFInTissueTCP, nFParTCP, nFTreatmentTCP, nFInTissueDim, nFInTum,
    //    nFInVes);
    //var1ParRange(kp, L, nFRefParInt, nFInTissueDim, nFInTum, nFInVes)
    //varArtTissue(P, nFArt, nFInTissueDim, nFRefParMean);
    //varArtTissue(N, P, nFDensInt, nFInTissueDim, nFRefParMean);
    //varArtTissue(N, P, nFDensInt, nFSigmaTumInt, nFInTissueDim, nFRefParMean);
    /*varErr(nFVarPar, nFMostRelPar, nFLeastRelPar, nFInTissueDim, nFInTum,
           nFInVes, L, P);*/
    //varParFromFiles(nFPar, nFInTissueDim, nFInTum, nFInVes);
    //varParFromFiles(nFInTissuePar, nFPar, nFTreatment);
    //varStoch(N, P, nFRefParInt, nFInTissueDim, nFInTum, nFInVes);
    writeInFiles(nFInTissuePar);
}



