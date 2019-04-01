#include "alloc.hpp"


/*------------------------------------------------------------------------------
 * This function allocates a 2D matrix nRow x nCol .
 *
 * Inputs:
 *  - nRow: number of rows of the tissue,
 *  - nCol: number of columns of the tissue.
 *
 * Outputs:
 *  - X: 2D matrix nRow x nCol.
------------------------------------------------------------------------------*/

double **alloc2D(const int nRow, const int nCol){
    double **X;
    X = new double *[nRow];
    for(int i(0); i < nRow; i++){
        X[i] = new double [nCol];
    }
    return X;
}


/*------------------------------------------------------------------------------
 * This function allocates a 3D matrix nLayer x nRow x nCol .
 *
 * Inputs:
 *  - nLayer: number of layers of the tissue,
 *  - nRow: number of rows of the tissue,
 *  - nCol: number of columns of the tissue.
 *
 * Outputs:
 *  - X: 3D matrix nLayer x nRow x nCol.
------------------------------------------------------------------------------*/

double ***alloc3D(const int nLayer, const int nRow, const int nCol){
    double ***X;
    X = new double **[nLayer];
    for(int i(0); i < nLayer; i++){
        X[i] = alloc2D(nRow, nCol);
    }
    return X;
}


/*------------------------------------------------------------------------------
 * This function deallocates a 2D matrix nRow x *.
 *
 * Inputs:
 *  - X: 2D matrix nRow x *,
 *  - nRow: number of rows of the tissue.
------------------------------------------------------------------------------*/

void free2D(double **X, const int nRow){
    for(int i(0); i < nRow; i++){
        delete [] X[i];
    }
    delete [] X;
}


/*------------------------------------------------------------------------------
 * This function deallocates a 3D matrix nLayer x nRow x *.
 *
 * Inputs:
 *  - X: 3D matrix nLayer x nRow x *,
 *  - nLayer: number of layers of the tissue,
 *  - nRow: number of rows of the tissue.
------------------------------------------------------------------------------*/

void free3D(double ***X, const int nLayer, const int nRow){
    for(int i(0); i < nLayer; i++){
        free2D(X[i], nRow);
    }
    delete [] X;
}
