#ifndef DEF_ALLOC
#define DEF_ALLOC

double **alloc2D(const int nRow, const int nCol);
double ***alloc3D(const int nLayer, const int nRow, const int nCol);
void free2D(double **X, const int nRow);
void free3D(double ***X, const int nLayer, const int nRow);

#endif
