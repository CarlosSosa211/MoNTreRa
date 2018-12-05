#include "alloc.hpp"

double **alloc2D(const int nRow, const int nCol){
    double **X;
    X = new double *[nRow];
    for(int i(0); i < nRow; i++){
        X[i] = new double [nCol];
    }
    return X;
}


double ***alloc3D(const int nLayer, const int nRow, const int nCol){
    double ***X;
    X = new double **[nLayer];
    for(int i(0); i < nLayer; i++){
        X[i] = alloc2D(nRow, nCol);
    }
    return X;
}


void free2D(double **X, const int nRow){
    for(int i(0); i < nRow; i++){
        delete [] X[i];
    }
    delete [] X;
}


void free3D(double ***X, const int nLayer, const int nRow){
    for(int i(0); i < nLayer; i++){
        free2D(X[i], nRow);
    }
    delete [] X;
}
