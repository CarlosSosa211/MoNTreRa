#include "createAndReadFiles.hpp"

using namespace std;


/*------------------------------------------------------------------------------
 * This function creates the initial configuration of an artificial tissue in
 * terms of tumour and endothelial cells.
 *
 * Inputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - tumDens: initial tumour density,
 *  - sigmaTum: initial tumour sigma. A small value (close to 0), indicates that
 *  tumour cells form a cluster at the center of the tissue. For a big value
 *  (close to 1), they appear to follow a uniform distribution.
 *  - vascDens: initial vascular density,
 *  - sigmaVasc: initial vascular sigma. A small value (close to 0), indicates
 *  that endothelial cells are placed on a regular grid. For a big value
 *  (close to 1), they appear to follow a uniform distribution.
 *  - vascDens: initial vascular density.
 *
 * Outputs:
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double sigmaTum,
                  const double vascDens, const double sigmaVasc,
                  vector<bool> &inTum, vector<bool> &inVes){
    int ivd, ivd2;
    int mindim, mindim2, sqrtmin, tumToDist, vesToDist;
    int nrowNcol, nrowNcolNlayer;
    vector<int> div;
    vector<double> diff;

    nrowNcol = nrow * ncol;
    nrowNcolNlayer = nrowNcol * nlayer;
    tumToDist = tumDens * nrowNcolNlayer;

    if(vascDens){
        mindim = min(nrow, ncol);
        mindim2 = mindim * mindim;
        sqrtmin = sqrt(mindim);

        for(int l(1); l < sqrtmin; l+=2){
            if(!(nrow % l) && !(ncol % l)){
                div.push_back(l);
                diff.push_back(fabs(1.0 / (l * l) - vascDens));
                div.push_back(mindim / l);
                diff.push_back(fabs(double(l * l) / double(mindim2) -
                                    vascDens));
            }
        }

        ivd = div.at(min_element(diff.begin(), diff.end()) - diff.begin());
        ivd2 = ivd * ivd;
        vesToDist = min(nrowNcol / ivd2, nrowNcolNlayer - tumToDist);
    }

    else{
        vesToDist = 0;
    }

    int halfIvd, halfNrow, halfNcol, halfNlayer;
    int imHalfNrow2, lmHalfNlayer2;
    int lnrowNcol;
    int irr2, ishift;
    double R, RR;
    vector<simpCell> map(nrowNcolNlayer);

    halfIvd = 0.5 * ivd;
    halfNrow = 0.5 * nrow;
    halfNcol = 0.5 * ncol;
    halfNlayer = 0.5 * nlayer;
    R = sqrt(halfNrow * halfNrow + halfNcol * halfNcol + halfNlayer *
             halfNlayer);
    RR = 1.001 * sqrt(2.0 * halfIvd * halfIvd);

    int k(0);

    for(int l(0); l < nlayer; l++){
        lmHalfNlayer2 = (l - halfNlayer) * (l - halfNlayer);
        lnrowNcol = l * nrowNcol;
        for(int i(0); i < nrow; i++){
            imHalfNrow2 = (i - halfNrow) * (i - halfNrow);
            irr2 = (i % ivd - halfIvd) * (i % ivd - halfIvd);
            ishift = i / ivd * ncol / ivd;
            for(int j(0); j < ncol; j++){
                map.at(k).r = sqrt(lmHalfNlayer2 + imHalfNrow2 +
                                   (j - halfNcol) * (j - halfNcol)) / R;
                map.at(k).rr = lnrowNcol + ishift + j / ivd +
                        sqrt(irr2 + (j % ivd - halfIvd) * (j % ivd - halfIvd)) /
                        RR;
                map.at(k).k = k;
                map.at(k).tum = 0;
                map.at(k).ves = 0;
                k++;
            }
        }
    }
    k = 0;

    sort(map.begin(), map.end(), compR);

    int m;
    double n;
    default_random_engine gen;
    normal_distribution<double> distTum(0, sigmaTum);
    bool tooBigTumDens(tumDens > 0.9);
    bool tooSmallSigmaTum(sigmaTum < 0.15);

    if(!tooBigTumDens && !tooSmallSigmaTum){
        while(tumToDist > 0){
            n = distTum(gen);
            if(n >= 0.0 && n < 1.0){
                m = n * nrowNcolNlayer;
                if(!map.at(m).tum){
                    map.at(m).tum = 1;
                    tumToDist--;
                }
            }
        }
    }

    else if(tooBigTumDens){
        for(int k(0); k < map.size(); k++){
            map.at(k).tum = 1;
        }
        reverse(map.begin(), map.end());
        int tumToRemove(nrowNcolNlayer - tumToDist);
        while(tumToRemove > 0){
            n = distTum(gen);
            if(n >= 0.0 && n < 1.0){
                m = n * nrowNcolNlayer;
                if(map.at(m).tum){
                    map.at(m).tum = 0;
                    tumToRemove--;
                }
            }
        }
    }

    else{
        bool cond;
        while(tumToDist > 0){
            n = distTum(gen);
            if(n >= 0.0 && n < 1.0){
                m = n * nrowNcolNlayer;
                cond = true;
                while(cond){
                    if(!map.at(m).tum){
                        map.at(m).tum = 1;
                        tumToDist--;
                        cond = false;
                    }
                    else{
                        m++;
                    }
                    if(m == nrowNcolNlayer){
                        m = 0;
                    }
                }
            }
        }
    }

    sort(map.begin(), map.end(), compRR);

    int nivd2;
    bool availSpace, firstRound(true);
    normal_distribution<double> distVes(0, sigmaVasc);

    k = 0;
    while(vesToDist > 0){
        availSpace = false;
        for(int j(0); j < ivd2; j++){
            if(!map.at(k + j).ves && !map.at(k + j).tum){
                availSpace = true;
                break;
            }
        }

        if(availSpace){
            n = distVes(gen);
            while(n < 0.0 || n >= 1.0){
                n = distVes(gen);
            }
            nivd2 = n * ivd2;
            m = k + nivd2;
            while(map.at(m).ves || map.at(m).tum){
                nivd2++;
                if(nivd2 >= ivd2){
                    nivd2 = 0;
                }
                m = k + nivd2;
            }
            map.at(m).ves = 1;
            vesToDist--;
            availSpace = false;
        }

        if(firstRound){
            k += ivd2;
        }
        else{
            k = rand() % (nrowNcol / ivd2) * ivd2;
        }
        if(k >= nrowNcol - ivd2){
            firstRound = false;
        }
    }

    sort(map.begin(), map.end(), compK);
    for(int k(0); k < nrowNcol; k++){
        if(map.at(k).ves == 1){
            for(int l(1); l < nlayer; l++){
                map.at(l * nrowNcol + k).ves = 1;
                map.at(l * nrowNcol + k).tum = 0;
            }
        }
    }

    for(k = 0; k < nrowNcolNlayer; k++){
        inTum.push_back(map.at(k).tum);
        inVes.push_back(map.at(k).ves);
    }

    return 0;
}


/*------------------------------------------------------------------------------
 * This function creates the initial configuration of a tissue in terms of
 * endothelial cells. It is supposed to contain no tumour cells.
 *
 * Inputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - vascDens: initial vascular density,
 *  - sigmaVasc: initial vascular sigma. A small value (close to 0), indicates
 *  that endothelial cells are placed on a regular grid. For a big value
 *  (close to 1), they appear to follow a uniform distribution.
 *
 * Outputs:
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double vascDens, const double sigmaVasc,
                  vector<bool> &inVes){
    int ivd, ivd2;
    int mindim, mindim2, sqrtmin, tumToDist, vesToDist;
    int nrowNcol, nrowNcolNlayer;
    vector<int> div;
    vector<double> diff;

    nrowNcol = nrow * ncol;
    nrowNcolNlayer = nrowNcol * nlayer;

    if(vascDens){
        mindim = min(nrow, ncol);
        mindim2 = mindim * mindim;
        sqrtmin = sqrt(mindim);

        for(int l(1); l < sqrtmin; l+=2){
            if(!(nrow % l) && !(ncol % l)){
                div.push_back(l);
                diff.push_back(fabs(1.0 / (l * l) - vascDens));
                div.push_back(mindim / l);
                diff.push_back(fabs(double(l * l) / double(mindim2) -
                                    vascDens));
            }
        }

        ivd = div.at(min_element(diff.begin(), diff.end()) - diff.begin());
        ivd2 = ivd * ivd;
        vesToDist = min(nrowNcol / ivd2, nrowNcolNlayer - tumToDist);
    }

    else{
        vesToDist = 0;
    }

    int halfIvd;
    int lnrowNcol;
    int irr2, ishift;
    double R, RR;
    vector<simpCell> map(nrowNcolNlayer);

    halfIvd = 0.5 * ivd;
    RR = 1.001 * sqrt(2.0 * halfIvd * halfIvd);

    int k(0);

    for(int l(0); l < nlayer; l++){
        lnrowNcol = l * nrowNcol;
        for(int i(0); i < nrow; i++){
            irr2 = (i % ivd - halfIvd) * (i % ivd - halfIvd);
            ishift = i / ivd * ncol / ivd;
            for(int j(0); j < ncol; j++){
                map.at(k).rr = lnrowNcol + ishift + j / ivd +
                        sqrt(irr2 + (j % ivd - halfIvd) * (j % ivd - halfIvd)) /
                        RR;
                map.at(k).k = k;
                map.at(k).ves = 0;
                k++;
            }
        }
    }

    sort(map.begin(), map.end(), compRR);

    bool availSpace, firstRound(true);
    int m, nivd2;
    double n;
    default_random_engine gen;
    normal_distribution<double> distVes(0, sigmaVasc);

    k = 0;
    while(vesToDist > 0){
        availSpace = false;
        for(int j(0); j < ivd2; j++){
            if(!map.at(k + j).ves){
                availSpace = true;
                break;
            }
        }

        if(availSpace){
            n = distVes(gen);
            while(n < 0.0 || n >= 1.0){
                n = distVes(gen);
            }
            nivd2 = n * ivd2;
            m = k + nivd2;
            while(map.at(m).ves){
                nivd2++;
                if(nivd2 >= ivd2){
                    nivd2 = 0;
                }
                m = k + nivd2;
            }
            map.at(m).ves = 1;
            vesToDist--;
            availSpace = false;
        }

        if(firstRound){
            k += ivd2;
        }
        else{
            k = rand() % (nrowNcol / ivd2) * ivd2;
        }
        if(k >= nrowNcol - ivd2){
            firstRound = false;
        }
    }

    sort(map.begin(), map.end(), compK);
    for(int k(0); k < nrowNcol; k++){
        if(map.at(k).ves == 1){
            for(int l(1); l < nlayer; l++){
                map.at(l * nrowNcol + k).ves = 1;
            }
        }
    }

    for(k = 0; k < nrowNcolNlayer; k++){
        inVes.push_back(map.at(k).ves);
    }

    return 0;
}


/*------------------------------------------------------------------------------
 * This function reads the files containing the dimensions and the initial
 * configuration in terms of tumour and endothelial cell of a tissue.
 *
 * Inputs:
 *  - nFInTissueDim: name of the file containing the dimensions of a tissue,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration of a tissue,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration of a tissue.
 *
 * Outputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial endothelial cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/

void readInFiles(const string nFInTissueDim, const string nFInTum,
                 const string nFInVes, int &nrow, int &ncol, int &nlayer,
                 double &cellSize, vector<bool> &inTum, vector<bool> &inVes){
    ifstream fInTissueDim(nFInTissueDim.c_str());

    fInTissueDim >> nrow >> ncol >> nlayer;
    fInTissueDim >> cellSize;
    fInTissueDim.close();

    ifstream fInTum(nFInTum.c_str());
    bool temp;

    fInTum >> temp;
    while(!fInTum.eof()){
        inTum.push_back(temp);
        fInTum >> temp;
    }
    fInTum.close();

    ifstream fInVes(nFInVes.c_str());

    fInVes >> temp;
    while(!fInVes.eof()){
        inVes.push_back(temp);
        fInVes >> temp;
    }
    fInVes.close();
}


/*------------------------------------------------------------------------------
 * This function reads the files containing the dimensions and the initial
 * configuration in terms of endothelial cell of a tissue.
 *
 * Inputs:
 *  - nFInTissueDim: name of the file containing the dimensions of a tissue,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration of a tissue.
 *
 * Outputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/

void readInFiles(const string nFInTissueDim, const string nFInVes, int &nrow,
                 int &ncol, int &nlayer, double &cellSize, vector<bool> &inVes){
    ifstream fInTissueDim(nFInTissueDim.c_str());

    fInTissueDim >> nrow >> ncol >> nlayer;
    fInTissueDim >> cellSize;
    fInTissueDim.close();

    ifstream fInVes(nFInVes.c_str());
    bool temp;

    fInVes >> temp;
    while(!fInVes.eof()){
        inVes.push_back(temp);
        fInVes >> temp;
    }
    fInVes.close();
}


/*------------------------------------------------------------------------------
 * This function reads the files containing the dimensions and the initial
 * vascular configuration parameters of a tissue.
 *
 * Inputs:
 *  - nFInTissueOxy: name of the file containing the dimensions of a tissue and
 * its initial vascular configuration parameters,
 *
 * Outputs:
 *  - art: 0, if the tissue is a histological specimen or 1 if it is artificial
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - vascDens: initial vascular density,
 *  - sigmaVasc: initial vascular sigma. A small value (close to 0), indicates
 *  that endothelial cells are placed on a regular grid. For a big value
 *  (close to 1), they appear to follow a uniform distribution.
------------------------------------------------------------------------------*/

void readInFilesOxy(const string nFInTissueOxy, bool &art, int &nrow, int &ncol,
                    int &nlayer, double &cellSize, double &vascDens,
                    double &sigmaVasc){
    ifstream fInTissueOxy(nFInTissueOxy.c_str());

    fInTissueOxy >> art;
    if(art){
        fInTissueOxy >> nrow >> ncol >> nlayer;
        fInTissueOxy >> cellSize;
        fInTissueOxy >> vascDens >> sigmaVasc;
    }
    fInTissueOxy.close();
}


/*------------------------------------------------------------------------------
 * This function reads the files containing the dimensions and the initial
 * configuration parameters of a tissue and the treatments considered for a tcp
 * simulation.
 *
 * Inputs:
 *  - nFInTissueTCP: name of the file containing the dimensions of a tissue and
 * its initial configuration parameters,
 *  - nFInTreatmentTCP: vector with the names of the files containing the
 *  treatments.
 *
 * Outputs:
 *  - art: 0, if the tissue is a histological specimen or 1 if it is artificial
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - tumDens: initial tumour density,
 *  - sigmaTum: initial tumour sigma. A small value (close to 0), indicates that
 *  tumour cells form a cluster at the center of the tissue. For a big value
 *  (close to 1), they appear to follow a uniform distribution.
 *  - vascDens: initial vascular density,
 *  - sigmaVasc: initial vascular sigma. A small value (close to 0), indicates
 *  that endothelial cells are placed on a regular grid. For a big value
 *  (close to 1), they appear to follow a uniform distribution,
 *  - treatment: vector containing the treatments.
------------------------------------------------------------------------------*/

void readInFilesTCP(const string nFInTissueTCP,
                    const vector<string> nFTreatmentTCP, bool &art, int &nrow,
                    int &ncol, int &nlayer, double &cellSize, double &tumDens,
                    double &sigmaTum, double &vascDens, double &sigmaVasc,
                    vector<Treatment> &treatment){
    ifstream fInTissueTCP(nFInTissueTCP.c_str());

    fInTissueTCP >> art;
    if(art){
        fInTissueTCP >> nrow >> ncol >> nlayer;
        fInTissueTCP >> cellSize;
        fInTissueTCP >> tumDens >> sigmaTum >> vascDens >> sigmaVasc;
        fInTissueTCP.close();
    }

    double fraction, totalDose, interval;
    int schedule;
    for(int i(0); i < nFTreatmentTCP.size(); i++){
        ifstream fTreatmentTCP(nFTreatmentTCP[i].c_str());
        fTreatmentTCP >> fraction >> totalDose >> interval >> schedule;;
        treatment.push_back(Treatment(fraction, totalDose, interval, schedule));
        fTreatmentTCP.close();
    }
}
