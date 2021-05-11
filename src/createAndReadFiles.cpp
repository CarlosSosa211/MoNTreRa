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
 *  - radRatioTum: initial tumour radius ratio. A small value (close to 0),
 *  indicates that tumour cells form a segment. For a value equal to 1, they
 *  form a disk,
 *  - sigmaTum: initial tumour sigma. A small value (close to 0), indicates that
 *  tumour cells form a cluster at the center of the tissue. For a big value
 *  (close to 1), they appear to follow a uniform distribution,
 *  - vascDens: initial vascular density,
 *  - sigmaVasc: initial vascular sigma. A small value (close to 0), indicates
 *  that endothelial cells are placed on a regular grid. For a big value
 *  (close to 1), they appear to follow a uniform distribution.
 *
 * Outputs:
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double radRatioTum,
                  const double sigmaTum, const double vascDens,
                  const double sigmaVasc, vector<bool> &inTum,
                  vector<bool> &inVes){
    int ivd, ivd2, l2;
    int mindim, mindim2, sqrtmin, tumToDist, vesToDist;
    int nrowNcol, nrowNcolNlayer;
    const double p2(radRatioTum * radRatioTum);
    vector<int> div;
    vector<double> diff;

    nrowNcol = nrow * ncol;
    nrowNcolNlayer = nrowNcol * nlayer;
    tumToDist = tumDens * nrowNcolNlayer;

    if(vascDens){
        mindim = min(nrow, ncol);
        mindim2 = mindim * mindim;
        sqrtmin = sqrt(mindim) + 1;

        for(int l(1); l < sqrtmin; l++){
            if(!(nrow % l) && !(ncol % l)){
                l2 = l * l;
                div.push_back(l);
                diff.push_back(fabs(1.0 / l2 - vascDens));
                div.push_back(mindim / l);
                diff.push_back(fabs(double(l2) / double(mindim2) - vascDens));
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
                                   p2 * (j - halfNcol) * (j - halfNcol)) / R;
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
        inTum.at(k) = map.at(k).tum;
        inVes.at(k) = map.at(k).ves;
    }

    return 0;
}


/*------------------------------------------------------------------------------
 * This function creates the initial configuration of an artificial tissue in
 * terms of tumour and endothelial cells, supposing a random uniform
 * distribution of the latter.
 *
 * Inputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - tumDens: initial tumour density,
 *  - radRatioTum: initial tumour radius ratio. A small value (close to 0),
 *  indicates that tumour cells form a segment. For a value equal to 1, they
 *  form a disk,
 *  - sigmaTum: initial tumour sigma. A small value (close to 0), indicates that
 *  tumour cells form a cluster at the center of the tissue. For a big value
 *  (close to 1), they appear to follow a uniform distribution,
 *  - vascDens: initial vascular density.
 *
 * Outputs:
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double radRatioTum,
                  const double sigmaTum, const double vascDens,
                  vector<bool> &inTum, vector<bool> &inVes){
    const int nrowNcol(nrow * ncol), nrowNcolNlayer(nrowNcol * nlayer);
    const double p2(radRatioTum * radRatioTum);
    int m, tumToDist, vesToDist;
    double n;
    vector<simpCell> map(nrowNcolNlayer);

    vesToDist = vascDens * nrowNcolNlayer;
    tumToDist = min(int(tumDens * nrowNcolNlayer), nrowNcolNlayer - vesToDist);

    int halfNrow, halfNcol, halfNlayer;
    int imHalfNrow2, lmHalfNlayer2;
    int lnrowNcol;
    double R;

    halfNrow = 0.5 * nrow;
    halfNcol = 0.5 * ncol;
    halfNlayer = 0.5 * nlayer;
    R = sqrt(halfNrow * halfNrow + halfNcol * halfNcol + halfNlayer *
             halfNlayer);

    int k(0);
    for(int l(0); l < nlayer; l++){
        lmHalfNlayer2 = (l - halfNlayer) * (l - halfNlayer);
        lnrowNcol = l * nrowNcol;
        for(int i(0); i < nrow; i++){
            imHalfNrow2 = (i - halfNrow) * (i - halfNrow);
            for(int j(0); j < ncol; j++){
                map.at(k).r = sqrt(lmHalfNlayer2 + imHalfNrow2 +
                                   p2 * (j - halfNcol) * (j - halfNcol)) / R;
                map.at(k).k = k;
                map.at(k).tum = 0;
                map.at(k).ves = 0;
                k++;
            }
        }
    }

    sort(map.begin(), map.end(), compR);

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

    sort(map.begin(), map.end(), compK);

    while(vesToDist > 0){
        n = double(rand()) / (RAND_MAX + 1.0);
        m = n * nrowNcolNlayer;
        if(!map.at(m).ves && !map.at(m).tum){
            map.at(m).ves = 1;
            vesToDist--;
        }
    }

    for(k = 0; k < nrowNcolNlayer; k++){
        inTum.at(k) = map.at(k).tum;
        inVes.at(k) = map.at(k).ves;
    }

    return 0;
}


/*------------------------------------------------------------------------------
 * This function creates the initial configuration of a tissue in terms of
 * endothelial cells, giving a tumour cell configuration.
 *
 * Inputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - vascDens: initial vascular density,
 *  - inTum: vector containing the initial tumour cell configuration.
 *
 * Outputs:
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double vascDens, const vector<bool> &inTum,
                  vector<bool> &inVes){
    const int nrowNcolNlayer(nrow * ncol * nlayer);
    vector<simpCell> map(nrowNcolNlayer);

    for(int k(0); k < nrowNcolNlayer; k++){
        map.at(k).tum = inTum.at(k);
        map.at(k).ves = 0;
    }

    int vesToDist(vascDens * nrowNcolNlayer);
    int m;
    double n;
    while(vesToDist > 0){
        n = double(rand()) / (RAND_MAX + 1.0);
        m = n * nrowNcolNlayer;
        if(!map.at(m).ves && !map.at(m).tum){
            map.at(m).ves = 1;
            vesToDist--;
        }
    }

    for(int k(0); k < nrowNcolNlayer; k++){
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
    int ivd, ivd2, l2;
    int mindim, mindim2, sqrtmin, vesToDist;
    int nrowNcol, nrowNcolNlayer;
    vector<int> div;
    vector<double> diff;

    nrowNcol = nrow * ncol;
    nrowNcolNlayer = nrowNcol * nlayer;

    if(vascDens){
        mindim = min(nrow, ncol);
        mindim2 = mindim * mindim;
        sqrtmin = sqrt(mindim) + 1;

        for(int l(1); l < sqrtmin; l++){
            l2 = l * l;
            if(!(nrow % l) && !(ncol % l)){
                div.push_back(l);
                diff.push_back(fabs(1.0 / (l2) - vascDens));
                div.push_back(mindim / l);
                diff.push_back(fabs(double(l2) / double(mindim2) - vascDens));
            }
        }

        ivd = div.at(min_element(diff.begin(), diff.end()) - diff.begin());
        ivd2 = ivd * ivd;
        vesToDist = nrowNcol / ivd2;
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
 * This function creates the initial configuration of an artificial tissue in
 * terms of tumour and endothelial cells, supposing a random uniform
 * distribution.
 *
 * Inputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - tumDens: initial tumour density,
 *  - vascDens: initial vascular density.
 *
 * Outputs:
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/

int createInFiles(const int nrow, const int ncol, const int nlayer,
                  const double tumDens, const double vascDens,
                  vector<bool> &inTum, vector<bool> &inVes){
    const int nrowNcolNlayer(nrow * ncol * nlayer);
    int m, tumToDist, vesToDist;
    double n;
    vector<simpCell> map(nrowNcolNlayer);
    vesToDist = vascDens * nrowNcolNlayer;
    tumToDist = min(int(tumDens * nrowNcolNlayer), nrowNcolNlayer - vesToDist);

    for(int k(0); k < nrowNcolNlayer; k++){
        map.at(k).tum = 0;
        map.at(k).ves = 0;
    }

    while(vesToDist > 0){
        n = double(rand()) / (RAND_MAX + 1.0);
        m = n * nrowNcolNlayer;
        if(!map.at(m).ves){
            map.at(m).ves = 1;
            vesToDist--;
        }
    }

    while(tumToDist > 0){
        n = double(rand()) / (RAND_MAX + 1.0);
        m = n * nrowNcolNlayer;
        if(!map.at(m).ves && !map.at(m).tum){
            map.at(m).tum = 1;
            tumToDist--;
        }
    }

    for(int k(0); k < nrowNcolNlayer; k++){
        inTum.at(k) = map.at(k).tum;
        inVes.at(k) = map.at(k).ves;
    }

    return 0;
}

/*------------------------------------------------------------------------------
 * This function creates the initial configuration of an artificial tissue in
 * terms of tumour and endothelial cells, supposing a random uniform
 * distribution for the latter.
 *
 * Inputs:
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - tumArea: initial tumour area,
 *  - radRatioTum: initial tumour radius ratio. A small value (close to 0),
 *  indicates that tumour cells form a segment. For a value equal to 1, they
 *  form a disk,
 *  - tumDens: initial tumour density in the vascular area,
 *  - vascDens: initial vascular density.
 *
 * Outputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - inTum: vector containing the initial tumour cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration.
------------------------------------------------------------------------------*/
int createInFiles(const double cellSize, const double tumArea,
                  const double radRatioTum, const double tumDens,
                  const double vascDens, int &nrow, int &ncol, int &nlayer,
                  vector<bool> &inTum, vector<bool> &inVes){
    const int r(sqrt(tumArea * radRatioTum / M_PI) / cellSize);
    int m, tumToDist, vesToDist;
    double n;

    nrow = 2.5 * r;
    ncol = 2.5 * r / radRatioTum;
    nlayer = 1;

    const int nrowNcol(nrow * ncol), nrowNcolNlayer(nrowNcol * nlayer);
    const double p2(radRatioTum * radRatioTum);
    int halfNrow, halfNcol, halfNlayer;
    int imHalfNrow2, lmHalfNlayer2;
    int lnrowNcol;
    double R;
    vector<simpCell> map(nrowNcolNlayer);

    halfNrow = 0.5 * nrow;
    halfNcol = 0.5 * ncol;
    halfNlayer = 0.5 * nlayer;
    R = sqrt(halfNrow * halfNrow + halfNcol * halfNcol + halfNlayer *
             halfNlayer);

    int k(0);
    for(int l(0); l < nlayer; l++){
        lmHalfNlayer2 = (l - halfNlayer) * (l - halfNlayer);
        lnrowNcol = l * nrowNcol;
        for(int i(0); i < nrow; i++){
            imHalfNrow2 = (i - halfNrow) * (i - halfNrow);
            for(int j(0); j < ncol; j++){
                map.at(k).r = sqrt(lmHalfNlayer2 + imHalfNrow2 +
                                   p2 * (j - halfNcol) * (j - halfNcol)) / R;
                map.at(k).k = k;
                map.at(k).tum = 0;
                map.at(k).ves = 0;
                k++;
            }
        }
    }

    sort(map.begin(), map.end(), compR);

    int nCellsTumArea(tumArea / (cellSize * cellSize));
    vesToDist = vascDens * nrowNcolNlayer;
    tumToDist = min(int(nCellsTumArea * tumDens),
                    nrowNcolNlayer - vesToDist);

    default_random_engine gen;
    normal_distribution<double> distTum(0, 0.33);

    while(tumToDist > 0){
        n = distTum(gen);
        if(n >= 0.0 && n < 1.0){
            m = n * nCellsTumArea;
            if(!map.at(m).tum){
                map.at(m).tum = 1;
                tumToDist--;
            }
        }
    }

    //    sort(map.begin(), map.end(), compK);
    //    reverse(map.begin(),map.end());

    //    std::normal_distribution<double> distVes(0, 1.0);

    while(vesToDist > 0){
        n = double(rand()) / (RAND_MAX + 1.0);
        m = n * nrowNcolNlayer;
        if(!map.at(m).ves && !map.at(m).tum){
            map.at(m).ves = 1;
            vesToDist--;
        }
    }

    //    while(vesToDist > 0){
    //        n = double(rand()) / (RAND_MAX + 1.0);
    //        m = n * nrowNcolNlayer;
    //        if(!map.at(m).ves){
    //            map.at(m).ves = 1;
    //            map.at(m).tum = 0;
    //            vesToDist--;
    //        }
    //    }

    //    while(vesToDist > 0){
    //        n = distVes(gen);
    //        if(n >= 0.0 && n < 1.0){
    //            m = n * nrowNcolNlayer;
    //            if(!map.at(m).ves){
    //                map.at(m).ves = 1;
    //                map.at(m).tum = 0;
    //                vesToDist--;
    //            }
    //        }
    //    }

    std::sort(map.begin(), map.end(), compK);

    inTum.resize(nrowNcolNlayer);
    inVes.resize(nrowNcolNlayer);

    for(k = 0; k < nrowNcolNlayer; k++){
        inTum.at(k) = map.at(k).tum;
        inVes.at(k) = map.at(k).ves;
    }

    return 0;
}

/*------------------------------------------------------------------------------
 * This function reads the files containing the dimensions and the initial
 * configuration in terms of tumour and endothelial cells of a tissue.
 *
 * Inputs:
 *  - nFInTissueDim: name of the file containing the dimensions of a tissue,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration of a tissue,
 *  - nFInVes: name of the file containing the initial endothelial cell
 *  configuration of a tissue,
 *  - nFInTreatment: name of the file containing the treatment.
 *
 * Outputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTum: vector containing the initial endothelial cell configuration,
 *  - inVes: vector containing the initial endothelial cell configuration,
 *  - treatment: treatment.
------------------------------------------------------------------------------*/

void readInFiles(const string nFInTissueDim, const string nFInTum,
                 const string nFInVes, int &nrow, int &ncol, int &nlayer,
                 double &cellSize, vector<bool> &inTum, vector<bool> &inVes,
                 const string nFTreatment, Treatment &treatment){
    ifstream fInTissueDim(nFInTissueDim.c_str());

    fInTissueDim >> nrow >> ncol >> nlayer >> cellSize;
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

    if(!nFTreatment.empty()){
        double fraction, totalDose, interval;
        int schedule;
        ifstream fTreatment(nFTreatment.c_str());
        fTreatment >> fraction >> totalDose >> interval >> schedule;
        treatment = Treatment(fraction, totalDose, interval, schedule);
    }
}


/*------------------------------------------------------------------------------
 * This function reads the files containing the dimensions and the initial
 * configuration in terms of tumour/endothelial cells of a tissue.
 *
 * Inputs:
 *  - nFInTissueDim: name of the file containing the dimensions of a tissue,
 *  - nFInTumVes: name of the file containing the initial tumour/endothelial
 *  cell configuration of a tissue.
 *
 * Outputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - inTumVes: vector containing the initial tumour/endothelial cell
 *  configuration.
------------------------------------------------------------------------------*/

void readInFiles(const string nFInTissueDim, const string nFInTumVes, int &nrow,
                 int &ncol, int &nlayer, double &cellSize,
                 vector<bool> &inTumVes){
    ifstream fInTissueDim(nFInTissueDim.c_str());

    fInTissueDim >> nrow >> ncol >> nlayer >> cellSize;
    fInTissueDim.close();

    ifstream fInTumVes(nFInTumVes.c_str());
    bool temp;

    fInTumVes >> temp;
    while(!fInTumVes.eof()){
        inTumVes.push_back(temp);
        fInTumVes >> temp;
    }
    fInTumVes.close();
}


/*------------------------------------------------------------------------------
 * This function reads the files containing the dimensions of a tissue.
 *
 * Inputs:
 *  - nFInTissueDim: name of the file containing the dimensions of a tissue,
 *
 * Outputs:
 *  - nrow: number of rows of the tissue,
 *  - ncol: number of columns of the tissue,
 *  - nlayer: number of layers of the tissue,
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue.
------------------------------------------------------------------------------*/

void readInFiles(const string nFInTissueDim, int &nrow, int &ncol, int &nlayer,
                 double &cellSize){
    ifstream fInTissueDim(nFInTissueDim.c_str());

    fInTissueDim >> nrow >> ncol >> nlayer;
    fInTissueDim >> cellSize;
    fInTissueDim.close();
}


/*------------------------------------------------------------------------------
 * This function reads the files containing the initial configuration
 * parameters of a tissue and a treatment.
 *
 * Inputs:
 *  - nFInTissuePar: name of the file containing the initial configuration
 *  parameters of a tissue,
 *  - nFInTreatment: name of the file containing the treatment.
 *
 * Outputs:
 *  - cellSize: length of the side of square cells, corresponding to a voxel
 *  of the tissue,
 *  - tumArea: initial tumour area,
 *  - radRatioTum: initial tumour radius ratio. A small value (close to 0),
 *  indicates that tumour cells form a segment. For a value equal to 1, they
 *  form a disk,
 *  - tumDens: initial tumour density in the vascular area,
 *  - vascDens: initial vascular density,
 *  - treatment: treatment.
------------------------------------------------------------------------------*/

void readInFiles(const string nFInTissuePar, double &cellSize, double &tumArea,
                 double &radRatioTum, double &tumDens, double &vascDens,
                 const string nFTreatment, Treatment &treatment){

    ifstream fInTissuePar(nFInTissuePar.c_str());
    fInTissuePar >> cellSize >> tumArea >> radRatioTum >> tumDens >> vascDens;
    fInTissuePar.close();

    if(!nFTreatment.empty()){
        double fraction, totalDose, interval;
        int schedule;
        ifstream fTreatmentTCP(nFTreatment.c_str());
        fTreatmentTCP >> fraction >> totalDose >> interval >> schedule;
        treatment = Treatment(fraction, totalDose, interval, schedule);
    }
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
 *  - radRatioTum: initial tumour radius ratio. A small value (close to 0),
 *  indicates that tumour cells form a segment. For a value equal to 1, they
 *  form a disk,
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
                    double &radRatioTum, double &sigmaTum, double &vascDens,
                    double &sigmaVasc, vector<Treatment> &treatment){
    ifstream fInTissueTCP(nFInTissueTCP.c_str());

    fInTissueTCP >> art;
    if(art){
        fInTissueTCP >> nrow >> ncol >> nlayer;
        fInTissueTCP >> cellSize;
        fInTissueTCP >> tumDens >> radRatioTum >> sigmaTum >> vascDens >>
                sigmaVasc;
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


/*------------------------------------------------------------------------------
 * This function reads a file containing the initial configuration parameters of
 * a tissue and creates and writes the files containing the corresponding
 * initial configurations of tumour and endothelial cells.
 *
 * Inputs:
 *  - nFInTissuePar: name of the file containing the initial configuration
 *  parameters of a tissue.
------------------------------------------------------------------------------*/

void writeInFiles(const string nFInTissuePar){
    int nrow, ncol, nlayer;
    double cellSize, tumArea, radRatioTum, tumDens, vascDens;
    vector<bool> inTum, inVes;

    readInFiles(nFInTissuePar, cellSize, tumArea, radRatioTum, tumDens,
                vascDens);
    createInFiles(cellSize, tumArea, radRatioTum, tumDens, vascDens, nrow, ncol,
                  nlayer, inTum, inVes);

    ofstream fInTissueDim("../InputFiles/tissueDim.dat");

    fInTissueDim << nrow << endl;
    fInTissueDim << ncol << endl;
    fInTissueDim << nlayer << endl;
    fInTissueDim << cellSize << endl;
    fInTissueDim.close();

    int nrowNcolNlayer(nrow * ncol * nlayer);
    ofstream fInTum("../InputFiles/inTum.dat");
    ofstream fInVes("../InputFiles/inVes.dat");

    for(int k(0); k < nrowNcolNlayer; k++){
        fInTum << inTum.at(k) << endl;
        fInVes << inVes.at(k) << endl;
    }

    fInTum.close();
    fInVes.close();
}


/*------------------------------------------------------------------------------
 * This function reads the files containing the dimensions and the initial
 * configuration of tumour cells of a tissue and creates and writes the file
 * containing the initial configuration of endothelial cells.
 *
 * Inputs:
 *  - nFInTissueDim: name of the file containing the dimensions of a tissue,
 *  - nFInTum: name of the file containing the initial tumour cell
 *  configuration of a tissue.
------------------------------------------------------------------------------*/

void writeInVesFile(const string nFInTissueDim, const std::string nFInTum){
    int nrow, ncol, nlayer;
    double cellSize, vascDens(0.038);
    vector<bool> inTum, inVes;

    readInFiles(nFInTissueDim, nFInTum, nrow, ncol, nlayer, cellSize, inTum);
    createInFiles(nrow, ncol, nlayer, vascDens, inTum, inVes);

    int nrowNcolNlayer(nrow * ncol * nlayer);
    ofstream fInVes("../InputFiles/inVes.dat");

    for(int k(0); k < nrowNcolNlayer; k++){
        fInVes << inVes.at(k) << endl;
    }
    fInVes.close();
}
