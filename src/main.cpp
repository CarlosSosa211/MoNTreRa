#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <time.h>
#include <vector>

#include "oxyTissue.hpp"
#include "rootSimulator.hpp"
#include "simpcell.hpp"
#include "tissue.hpp"
#include "treatment.hpp"

using namespace std;

double **alloc2D(const int nRow, const int nCol);
double ***alloc3D(const int nLayer, const int nRow, const int nCol);
int createInFiles(int nrow, int ncol, int nlayer,
                  double tumDens, double sigmaTum,
                  double vascDens, double sigmaVasc,
                  vector<bool> &inTum, vector<bool> &inVes);
void evalR(const int nMethod, const int nModel);
void free2D(double **X, const int nRow);
void free3D(double ***X, const int nLayer, const int nRow);
void model(double *x, double *y);
void morris(const int K, const int L, const int N, const int nOut,
	    const int p, const double *x0, const double *h,
	    double **mu, double **sigma);
void morrisRT(const int N, const int p, const string nFRefParInt);
void morrisToy(const int N, const int p, const string nFRefParInt);
void morrisVarRange(const int K, const int kp, const int L,
		    const int N, const int nOut, const int p,
		    const string nFRefParInt,
		    double ***mu, double ***sigma);
void morrisVarRangeRT(const int kp, const int L, const int N, const int p,
		      const string nFRefParInt);
void morrisVarRangeToy(const int kp, const int L, const int N, const int p,
		       const string nFRefParInt);
void sobol(const int K, const int N, const int nOut,
	   const double *x0, const double *h,
	   double **SI, double **TSI,
	   double ***SIConv, double ***TSIConv);
void sobolFromFiles(int K);
void sobolRT(const int N, const string nFRefParInt);
void sobolToy(const int N, const string nFRefParInt);
void toyModel(double *x, double *y);


int main(){
  //const int N(1e5);
  const int kp(0), L(10), p(20), N(100);
  //const int nMethod(1), nModel(0);
  string nFRefParInt("../InputFiles/refParIntRT.dat");
  //string nFRefParInt("../InputFiles/refParIntToy.dat");

  srand(time(NULL));
  //evalR(nMethod, nModel);
  morrisRT(N, p, nFRefParInt);
  //morrisToy(N, p, nFRefParInt);
  //morrisVarRangeRT(kp, L, N, p, nFRefParInt);
  //morrisVarRangeToy(kp, L, N, p, nFRefParInt);
  //sobolRT(N, nFRefParInt);
  //sobolToy(N, nFRefParInt);
  //sobolFromFiles(2);

  return 0;
}


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


int createInFiles(int nrow, int ncol, int nlayer,
                  double tumDens, double sigmaTum,
                  double vascDens, double sigmaVasc,
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
	diff.push_back(fabs(double(l * l) / double(mindim2) - vascDens));
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
  R = sqrt(halfNrow * halfNrow + halfNcol * halfNcol + halfNlayer * halfNlayer);
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
	  sqrt(irr2 + (j % ivd - halfIvd) * (j % ivd - halfIvd)) / RR;
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
    const int nOut(4);
    double y[nOut];  
    ofstream fYTumDens("../OutputFiles/YTumDens.res");
    ofstream fYIntTumDens("../OutputFiles/YIntTumDens.res");
    ofstream fYTimeTo95("../OutputFiles/YTimeTo95.res");
    ofstream fYTimeTo99("../OutputFiles/YTimeTo99.res");
    
    for(int i(0); i < nEv; i++){
      for(int k(0); k < K; k++){
	fX >> x[k];
      }
      model(x, y);
      //cout << i + 1 << " out of " << nEv << " evaluations of the model" << endl;
      //cout << "---------------------------------------------" << endl;

      fYTumDens    << y[0] << endl;
      fYIntTumDens << y[1] << endl;
      fYTimeTo95   << y[2] << endl;
      fYTimeTo99   << y[3] << endl;
    }

    fYTumDens.close();
    fYIntTumDens.close();
    fYTimeTo95.close();
    fYTimeTo99.close();
    break;
  }
  }

  fX.close();
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

  
void model(double *x, double *y){
  int k(0);
  /*int nrow(90), ncol(90), nlayer(1);
    cout << "Creating tissue with: "  << endl;
    cout << nrow << " row";
    if(nrow > 1){
    cout << "s";
    }
    cout << endl;
    cout << ncol << " column";
    if(ncol > 1){
    cout << "s";
    }
    cout << endl;
    cout << nlayer << " layer";
    if(nlayer > 1){
    cout << "s";
    }
    cout << endl;

    cout << "---------------------------------------------" << endl;*/

  /*double tumDensRef(0.4), sigmaTumRef(0.1);
    double vascDensRef(0.03), sigmaVascRef(0.1);

    double tumDens   = tumDensRef;
    double sigmaTum  = sigmaTumRef;
    double vascDens  = vascDensRef;
    double sigmaVasc = sigmaVascRef;

    cout << "tumDens: " << tumDens << endl;
    cout << "sigmaTum: " << sigmaTum << endl;
    cout << "vascDens: " << vascDens << endl;
    cout << "sigmaVasc: " << sigmaVasc << endl;*/

  int nrow, ncol, nlayer;
  double cellSize;

  std::ifstream fInTissueDim("../InputFiles/tissueDim.dat");

  fInTissueDim >> nrow >> ncol >> nlayer;
  fInTissueDim >> cellSize;
  fInTissueDim.close();

  vector<bool> inTum;
  std::ifstream fInTum("../InputFiles/inTum.dat");
  bool temp;

  fInTum >> temp;
  while(!fInTum.eof()){
    inTum.push_back(temp);
    fInTum >> temp;
  }

  fInTum.close();

  vector<bool> inVes;
  std::ifstream fInVes("../InputFiles/inVes.dat");

  fInVes >> temp;
  while(!fInVes.eof()){
    inVes.push_back(temp);
    fInVes >> temp;
  }
  fInVes.close();

  /*createInFiles(nrow, ncol, nlayer, tumDens, sigmaTum,
    vascDens, sigmaVasc, inTum, inVes);*/

  vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
  const double tumGrowth(1.0);
  const int edgeOrder(1);
  vector<double> cycDur = {0.55, 0.2, 0.15, 0.1};
  const double res(1.0);
  const double ang(1.0);

  const double tumTime(x[k]);
  k++;
  const double fibTime(x[k]);
  k++;
  const double vascTumTime(x[k]);
  k++;
  const double vegfThres(x[k]);
  k++;
  vector<double> alpha(8), beta(8);
  for(int i(0); i < 8; i++){
    alpha[i] = x[k];
    k++;
    beta[i]  = alpha[i] / x[k];
    k++;
  }
  const double doseThres(x[k]);
  k++;
  const double arrestTime(x[k]);
  k++;
  const double hypNecThres(x[k]);
  k++;
  const double dose(x[k]);
  k++;
  double DO2(x[k]);
  k++;
  double Vmax(x[k]);
  k++;
  const double Km(x[k]);
  k++;
  double Dvegf(x[k]);
  k++;
  double VmaxVegf(x[k]);
  k++;
  const double KmVegf(x[k]);
  k++;
  const double pO2NormVes(x[k]);
  k++;
  const double pO2TumVes(x[k]);
  k++;
  const double hypThres(x[k]);
  k++;
  const double hypVegf(x[k]);
  k++;

  cout << "tumTime: " << tumTime << " h" <<endl;
  cout << "fibTime: " << fibTime << " h" << endl;
  cout << "vascTumTime: " << vascTumTime << " h" << endl;
  cout << "vegfThres: " << vegfThres << " mol/um^3" << endl;
  cout << "alpha : ";
  for(int i(0); i < 8; i++){
    cout << alpha[i] << " ";
  }
  cout << "Gy^-1" << endl;
  cout << "beta : ";
  for(int i(0); i < 8; i++){
    cout << beta[i] << " ";
  }
  cout << "Gy^-2" << endl;
  cout << "doseThres: " << doseThres << " Gy" << endl;
  cout << "arrestTime: " << arrestTime << " h" << endl;
  cout << "hypNecThres: " << hypNecThres <<" mmHg" << endl;
  cout << "dose: " << dose << " Gy" << endl;
  cout << "DO2: " << DO2 << " um^2/ms" << endl;
  cout << "Vmax: " << Vmax << " mmHg/ms" << endl;
  cout << "Km: " << Km << " mmHg" << endl;
  cout << "Dvegf: " << Dvegf << " um^2/ms" << endl;
  cout << "VmaxVegf: " << VmaxVegf << " mol/um^3ms" << endl;
  cout << "KmVegf: " << KmVegf << " mol/um^3" << endl;
  cout << "pO2NormVes: " << pO2NormVes << " mmHg" << endl;
  cout << "pO2TumVes: " << pO2TumVes << " mmHg" << endl;
  cout << "hypThres: " << hypThres << " mmHg" << endl;
  cout << "hypVegf: " << hypVegf << " mol/um^3" << endl;


  Treatment *treatment;
  treatment = new Treatment(dose, 80.0, 24.0, 0);

  Coupler *coupler;
  Tissue *model1;
  OxyTissue *model2;

  model1 = new Tissue(nrow, ncol, nlayer, cellSize, inTum, inVes,
		      tumGrowth, tumTime, edgeOrder, cycDur, cycDistrib,
		      res, fibTime, ang, vascTumTime, vegfThres, alpha,
		      beta, doseThres, arrestTime, treatment, hypNecThres);

  double simTimeStep, oxySimTimeStep, sclFac, simTime;

  simTimeStep    = 6.0; //h;
  oxySimTimeStep = 10.0; //ms;
  sclFac = 3.6e6 * simTimeStep / oxySimTimeStep;
    
  DO2 *= oxySimTimeStep;
  Vmax *= oxySimTimeStep;
  Dvegf *= oxySimTimeStep;
  VmaxVegf *= oxySimTimeStep;


  model2 = new OxyTissue(nrow, ncol, nlayer, cellSize, inVes,
			 Dvegf, DO2, Vmax, Km, pO2NormVes, pO2TumVes,
			 hypThres, VmaxVegf, KmVegf, hypVegf);

  coupler = new Coupler(model1, model2);

  RootSimulator *sim;

  sim = new RootSimulator(coupler, simTimeStep,
			  oxySimTimeStep, sclFac);
  simTime = treatment->getDuration();
  sim->initSim();
  sim->simulate(simTimeStep, simTime);
  sim->stop();

  double finTumDens(model1->getOut()->at(0));
  double timeTo95(model1->getOut()->at(12));
  double timeTo99(model1->getOut()->at(13));
  double intTumDens(model1->getOut()->at(22));

  delete treatment;
  delete model1;
  delete model2;
  delete coupler;
  delete sim;

  cout << "finTumDens: " << finTumDens << endl;
  cout << "intTumDens: " << intTumDens << endl;
  cout << "timeTo95: " << timeTo95 << " h" << endl;
  cout << "timeTo99: " << timeTo99 << " h" << endl;
  
  y[0] = finTumDens;
  y[1] = intTumDens;
  y[2] = timeTo95;
  y[3] = timeTo99;
}


void morris(const int K, const int L, const int N, const int nOut,
	    const int p, const double *x0, const double *h,
	    double **mu, double **sigma){
  vector<int> vP;

  for(int k(0); k < K; k++){
    vP.push_back(k);
  }

  int diffk, nEv(0);
  const int M(K + 1);
  const int nEvTot(L * (K + 1) * N);
  const double _pm1(1.0 / (p - 1.0));
  const double delta(0.5 * _pm1 * p);
  const double _delta(1.0 / delta);
  double xp[M];
  double **B, **Bp;
  double **y;
  double ***elEff;

  B     = alloc2D(M, K);
  Bp    = alloc2D(M, K);
  y     = alloc2D(M, nOut);
  elEff = alloc3D(nOut, K, N);
  
  bool perm;

  for(int n(0); n < N; n++){
    for(int m(0); m < M; m++){
      xp[m] = _pm1 * (rand() % ((p - 2) / 2 + 1));
      for(int k(0); k < K; k++){
	B[m][k] = 0.0;
      }
    }
    
    for(int m(1); m < M; m++){
      for(int k(0); k < m ; k++){
	B[m][k] = 1.0;
      }
    }
    
    for(int k(0); k < K; k++){    
      perm = rand() % 2;
      for(int m(0); m < M; m++){
	if(perm){
	  B[m][k] = 1.0 - B[m][k];
	}
	B[m][k] = xp[k] + delta * B[m][k];
      }
    }

    random_shuffle(vP.begin(), vP.end());

    for(int m(0); m < M; m++){
      for(int k(0); k < K; k++){
	Bp[m][k] = x0[k] + B[m][vP[k]] * h[k];
      }
    }

    for(int m(0); m < M; m++){
      //toyModel(Bp[m], y[m]);
      model(Bp[m], y[m]);
      nEv++;
      cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
      cout << "---------------------------------------------" << endl;
    }

    for(int m(1); m < M; m++){
      diffk = find(vP.begin(), vP.end(), m - 1) - vP.begin();
      for(int j(0); j < nOut; j++){
	elEff[j][diffk][n] = fabs(_delta * (y[m][j] - y[m - 1][j]));	
	mu[j][diffk] += elEff[j][diffk][n];
      }
    }
  }

  free2D(B, M);
  free2D(Bp, M);
  free2D(y, M);
  
  for(int k(0); k < K; k++){
    for(int j(0); j < nOut; j++){
      mu[j][k] /= N;
    }
  }

  for(int k(0); k < K; k++){
    for(int n(0); n < N; n++){
      for(int j(0); j < nOut; j++){
	sigma[j][k] += (elEff[j][k][n] - mu[j][k]) *
	  (elEff[j][k][n] - mu[j][k]);
      }
    }

    for(int j(0); j < nOut; j++){
      sigma[j][k] = sqrt(sigma[j][k] / (N - 1.0));
    }
  }

  free3D(elEff, nOut, K);
}


void morrisRT(const int N, const int p, const string nFRefParInt){
  const int K(34), nOut(4);
  double h[K], x0[K];
  ifstream fRefParInt(nFRefParInt.c_str());
      
  for(int k(0); k < K; k++){
    fRefParInt >> x0[k];
    fRefParInt >> h[k];
    h[k] -= x0[k];
  }
  
  fRefParInt.close();

  double **mu, **sigma;

  mu    = alloc2D(nOut, K);
  sigma = alloc2D(nOut, K);
  
  for(int j(0); j < nOut; j++){
    for(int k(0); k < K; k++){
      mu[j][k]    = 0.0;
      sigma[j][k] = 0.0;
    }
  }

  morris(K, 1, N, nOut, p, x0, h, mu, sigma);
  
  ofstream fMorrisTumDens("../OutputFiles/morrisTumDens.res");
  ofstream fMorrisIntTumDens("../OutputFiles/morrisIntTumDens.res");
  ofstream fMorrisTimeTo95("../OutputFiles/morrisTimeTo95.res");
  ofstream fMorrisTimeTo99("../OutputFiles/morrisTimeTo99.res");

  for(int k(0); k < K; k++){
    fMorrisTumDens    << mu[0][k] << " " << sigma[0][k] << endl;
    fMorrisIntTumDens << mu[1][k] << " " << sigma[1][k] << endl;
    fMorrisTimeTo95   << mu[2][k] << " " << sigma[2][k] << endl;
    fMorrisTimeTo99   << mu[3][k] << " " << sigma[3][k] << endl;
  }
  
  fMorrisTumDens.close();
  fMorrisIntTumDens.close();
  fMorrisTimeTo95.close();
  fMorrisTimeTo99.close();

  free2D(mu, nOut);
  free2D(sigma, nOut);
}


void morrisToy(const int p, const int N, const string nFRefParInt){
  const int K(5), nOut(1);
  double h[K], x0[K];
  ifstream fRefParInt(nFRefParInt.c_str());
      
  for(int k(0); k < K; k++){
    fRefParInt >> x0[k];
    fRefParInt >> h[k];
    h[k] -= x0[k];
  }
  
  fRefParInt.close();

  double **mu, **sigma;

  mu    = alloc2D(nOut, K);
  sigma = alloc2D(nOut, K);

  for(int j(0); j < nOut; j++){
    for(int k(0); k < K; k++){
      mu[j][k]    = 0.0;
      sigma[j][k] = 0.0;
    }
  }

  morris(K, 1, N, nOut, p, x0, h, mu, sigma);
  
  ofstream fMorrisY("../OutputFiles/morrisY.res");

  for(int k(0); k < K; k++){
    fMorrisY << mu[0][k] << " " << sigma[0][k] << endl;
  }
  
  fMorrisY.close();

  free2D(mu, nOut);
  free2D(sigma, nOut);
}


void morrisVarRange(const int K, const int kp, const int L,
		    const int N, const int nOut, const int p,
		    const string nFRefParInt,
		    double ***mu, double ***sigma){
  double h[K], x0[K];
  ifstream fRefParInt(nFRefParInt.c_str());
  
  for(int k(0); k < K; k++){
    fRefParInt >> x0[k];
    fRefParInt >> h[k];
    h[k] -= x0[k];
  }
  fRefParInt.close();

  int nEv(0);
  const int nEvTot((K + 1) * N * L);
  const double h10(0.1 * h[kp]);
  ofstream fVarRange("../OutputFiles/morrisVarRange.res");
  
  fVarRange << kp << endl;
  
  for(int l(0); l < L; l++){
    fVarRange << h[kp] << endl;
    morris(K, L, N, nOut, p, x0, h, mu[l], sigma[l]);
    h[kp] += h10;
  }

  fVarRange.close();
}


void morrisVarRangeRT(const int kp, const int L, const int N, const int p,
		      const string nFRefParInt){
  const int K(34), nOut(4);
  double ***mu, ***sigma;

  mu    = alloc3D(L, nOut, K);
  sigma = alloc3D(L, nOut, K);

  for(int l(0); l < L; l++){
    for(int j(0); j < nOut; j++){
      for(int k(0); k < K; k++){
	mu[l][j][k]    = 0.0;
	sigma[l][j][k] = 0.0;
      }
    }
  }

  morrisVarRange(K, kp, L, N, nOut, p, nFRefParInt, mu, sigma);

  for(int l(0); l < L; l++){  
    ofstream fMorrisTumDens("../OutputFiles/morrisTumDens_" + to_string(l) + ".res");
    ofstream fMorrisIntTumDens("../OutputFiles/morrisIntTumDens_" + to_string(l) + ".res");
    ofstream fMorrisTimeTo95("../OutputFiles/morrisTimeTo95_" + to_string(l) + ".res");
    ofstream fMorrisTimeTo99("../OutputFiles/morrisTimeTo99_" + to_string(l) + ".res");

    for(int k(0); k < K; k++){
      fMorrisTumDens    << mu[l][0][k] << " " << sigma[l][0][k] << endl;
      fMorrisIntTumDens << mu[l][1][k] << " " << sigma[l][1][k] << endl;
      fMorrisTimeTo95   << mu[l][2][k] << " " << sigma[l][2][k] << endl;
      fMorrisTimeTo99   << mu[l][3][k] << " " << sigma[l][3][k] << endl;
    }
    fMorrisTumDens.close();
    fMorrisIntTumDens.close();
    fMorrisTimeTo95.close();
    fMorrisTimeTo99.close();
  }

  free3D(mu, L, nOut);
  free3D(sigma, L, nOut);  
}


void morrisVarRangeToy(const int kp, const int L, const int N, const int p,
		       const string nFRefParInt){
  const int K(5), nOut(1);
  double ***mu, ***sigma;

  mu    = alloc3D(L, nOut, K);
  sigma = alloc3D(L, nOut, K);

  for(int l(0); l < L; l++){
    for(int j(0); j < nOut; j++){
      for(int k(0); k < K; k++){
	mu[l][j][k]    = 0.0;
	sigma[l][j][k] = 0.0;
      }
    }
  }

  morrisVarRange(K, kp, L, N, nOut, p, nFRefParInt, mu, sigma);

  for(int l(0); l < L; l++){  
    ofstream fMorrisY("../OutputFiles/morrisY_" + to_string(l) + ".res");
    for(int k(0); k < K; k++){
      fMorrisY << mu[l][0][k] << " " << sigma[l][0][k] << endl;
    }
    fMorrisY.close();
  }

  free3D(mu, L, nOut);
  free3D(sigma, L, nOut);  
}


void sobol(const int K, const int N, const int nOut,
	   const double *x0, const double *h,
	   double **SI, double **TSI,
	   double ***SIConv, double ***TSIConv){
  int nEv(0), nEvTot((K + 2) * N);
  double **Xa, **Xb, **Xc;
  double **Ya, **Yb, **Yc;

  Xa = alloc2D(N, K);
  Xb = alloc2D(N, K);
  Xc = alloc2D(N, K);

  Ya = alloc2D(N, nOut);
  Yb = alloc2D(N, nOut);
  Yc = alloc2D(N, nOut);

  for(int i(0); i < N; i++){
    //Matrices for a model considering the Legendre polynomial of degree d
    //Xa[i][0] = rand() % 5 + 1;
    //Xa[i][1] = -1.0 + 2.0 * (double)(rand()) / (double)(RAND_MAX);
    //Xb[i][0] = rand() % 5 + 1;
    //Xb[i][1] = -1.0 + 2.0 * (double)(rand()) / (double)(RAND_MAX);
    for(int k(0); k < K; k++){
      Xa[i][k] = x0[k] + (double)(rand()) / (double)(RAND_MAX) * h[k];
      Xb[i][k] = x0[k] + (double)(rand()) / (double)(RAND_MAX) * h[k];
    }
  }

  for(int i(0); i < N; i++){
    toyModel(Xa[i], Ya[i]);
    //model(Xa[i], Ya[i]);
    nEv++;
    //cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
    //cout << "---------------------------------------------" << endl;
    toyModel(Xb[i], Yb[i]);
    //model(Xb[i], Yb[i]);
    nEv++;
    //cout << nEv << " out of " << nEvTot << " evaluations of the model" << endl;
    //cout << "---------------------------------------------" << endl;
  }

  int iConv(0), nConv;
  double alpha[nOut], beta[nOut], sigma2[nOut], f0[nOut];
  double f02, f02Conv, sigma2Conv;
  
  for(int i(0); i < N; i++){
    for(int k(0); k < K; k++){
      Xc[i][k] = Xa[i][k];
    }
  }

  double _N(1.0 / N);

  nConv = 2;
  
  for(int j(0); j < nOut; j++){  
    alpha[j]  = 0.0;
    beta[j]   = 0.0;
    f0[j]     = 0.0;
    sigma2[j] = 0.0;
  }
  
  for(int i(0); i < N; i++){
    Xc[i][0] = Xb[i][0];
    
    toyModel(Xc[i], Yc[i]);
    //model(Xc[i], Yc[i]);
    nEv++;
    //cout << nEv << " out of " << nEvTot << " evaluations of the model";
    //cout << "---------------------------------------------" << endl;
    
    for(int j(0); j < nOut; j++){  
      alpha[j]  += Yb[i][j] * Yc[i][j];
      beta[j]   += (Ya[i][j] - Yc[i][j]) * (Ya[i][j] - Yc[i][j]);
      f0[j]     += Yb[i][j] + Yc[i][j];
      sigma2[j] += Yb[i][j] * Yb[i][j] + Yc[i][j] * Yc[i][j];
    }
	
    if(i == nConv){
      for(int j(0); j < nOut; j++){  
	f02Conv    = 0.25 / (nConv * nConv) * f0[j] * f0[j] ;
	sigma2Conv = 0.5 / nConv * sigma2[j] - f02Conv;
	SIConv[iConv][j][0]  = (alpha[j] / nConv - f02Conv) / sigma2Conv;
	TSIConv[iConv][j][0] = 0.5 * beta[j] / (nConv * sigma2Conv);
      }
      iConv++;
      nConv *= 2;
    }
  }
  
  for(int j(0); j < nOut; j++){
    f02 = 0.25 * _N * _N * f0[j] * f0[j];
    sigma2[j] = 0.5 * _N * sigma2[j] - f02;
    SI[j][0] = (_N * alpha[j] - f02) / sigma2[j];
    TSI[j][0] = 0.5 * _N *  beta[j] / sigma2[j];
  }

  for(int k(1); k < K; k++){
    iConv = 0;
    nConv = 2;
    
    for(int j(0); j < nOut; j++){  
      alpha[j]  = 0.0;
      beta[j]   = 0.0;
      f0[j]     = 0.0;
      sigma2[j] = 0.0;
    }
    
    for(int i(0); i < N; i++){
      Xc[i][k - 1] = Xa[i][k - 1];
      Xc[i][k] = Xb[i][k];
      
      toyModel(Xc[i], Yc[i]);
      //model(Xc[i], Yc[i]);
      nEv++;
      //cout << nEv << " out of " << nEvTot << " evaluations of the model";
      //cout << "---------------------------------------------" << endl;
      
      for(int j(0); j < nOut; j++){  
	alpha[j]  += Yb[i][j] * Yc[i][j];
	beta[j]   += (Ya[i][j] - Yc[i][j]) * (Ya[i][j] - Yc[i][j]);
	f0[j]     += Yb[i][j] + Yc[i][j];
	sigma2[j] += Yb[i][j] * Yb[i][j] + Yc[i][j] * Yc[i][j];
      }
	
      if(i == nConv){
	for(int j(0); j < nOut; j++){ 
	  f02Conv    = 0.25 / (nConv * nConv) * f0[j] * f0[j] ;
	  sigma2Conv = 0.5 / nConv * sigma2[j] - f02Conv;
	  SIConv[iConv][j][k]  = (alpha[j] / nConv - f02Conv) / sigma2Conv;
	  TSIConv[iConv][j][k] = 0.5 * beta[j] / (nConv * sigma2Conv);
	}
	iConv++;
	nConv *= 2;
      }
    }
    
    for(int j(0); j < nOut; j++){
      f02       = 0.25 * _N * _N * f0[j] * f0[j] ;
      sigma2[j] = 0.5 * _N * sigma2[j] - f02;
      SI[j][k]  = (_N * alpha[j] - f02) / sigma2[j];
      TSI[j][k] = 0.5 * _N *  beta[j] / sigma2[j];
    }
  }

  free2D(Xa, N);
  free2D(Xb, N);
  free2D(Xc, N);

  free2D(Ya, N);
  free2D(Yb, N);
  free2D(Yc, N);
}


void sobolFromFiles(int K){
  int N;
  vector<double> Ya, Yb, Yc;

  ifstream fYa("../OutputSobol/Ya.txt");
  double temp;
  fYa >> temp;
  while(!fYa.eof()){
    Ya.push_back(temp);
    fYa >> temp;
  }
  fYa.close();

  N = Ya.size();

  ifstream fYb("../OutputSobol/Yb.txt");
  fYb >> temp;
  while(!fYb.eof()){
    Yb.push_back(temp);
    fYb >> temp;
  }
  fYb.close();

  double alpha, beta, sigma2, f0;
  vector<double> SI(K), TSI(K);
  ofstream fSens("../OutputFiles/sobol.res");

  fSens << K << endl;

  for(int k(0); k < K; k++){
    ifstream fYc("../OutputSobol/Yc" + to_string(k) + ".txt");
    Yc.resize(N);
    alpha = 0.0;
    beta = 0.0;
    sigma2 = 0.0;
    f0 = 0.0;
    for(int i(0); i < N; i++){
      fYc >> Yc[i];
      alpha += Yb[i] * Yc[i];
      beta += (Ya[i] - Yc[i]) * (Ya[i] - Yc[i]);
      f0 += Yb[i] + Yc[i];
      sigma2 += Yb[i] * Yb[i] + Yc[i] * Yc[i];
    }
    alpha /= N;
    f0 /= 2.0 * N;
    sigma2 /= 2.0 * N;
    sigma2 -= f0 * f0;
    SI[k] = (alpha - f0 * f0) / sigma2;
    TSI[k] = beta / (2.0 * N * sigma2);
    fSens << SI[k] << " " << TSI[k] << endl;
    fYc.close();
  }

  fSens.close();
}

void sobolRT(const int N, const string nFRefParInt){
  const int K(34), NConv(log(N) / log(2.0)), nOut(4);
  double h[K], x0[K];
  ifstream fRefParInt(nFRefParInt.c_str());
  
  for(int k(0); k < K; k++){
    fRefParInt >> x0[k];
    fRefParInt >> h[k];
    h[k] -= x0[k];
  }
  
  fRefParInt.close();

  double **SI, **TSI;
  double ***SIConv, ***TSIConv;
  
  SI  = alloc2D(nOut, K);
  TSI = alloc2D(nOut, K);
  SIConv  = alloc3D(NConv, nOut, K);
  TSIConv = alloc3D(NConv, nOut, K);

  sobol(K, N, nOut, x0, h, SI, TSI, SIConv, TSIConv);

  ofstream fSobolTumDens("../OutputFiles/sobolTumDens.res");
  ofstream fSobolIntTumDens("../OutputFiles/sobolIntTumDens.res");
  ofstream fSobolTimeTo95("../OutputFiles/sobolTimeTo95.res");
  ofstream fSobolTimeTo99("../OutputFiles/sobolTimeTo99.res");
  ofstream fConvSITumDens("../OutputFiles/convSITumDens.res");
  ofstream fConvTSITumDens("../OutputFiles/convTSITumDens.res");
  ofstream fConvSIIntTumDens("../OutputFiles/convSIIntTumDens.res");
  ofstream fConvTSIIntTumDens("../OutputFiles/convTSIIntTumDens.res");
  ofstream fConvSITimeTo95("../OutputFiles/convSITimeTo95.res");
  ofstream fConvTSITimeTo95("../OutputFiles/convTSITimeTo95.res");
  ofstream fConvSITimeTo99("../OutputFiles/convSITimeTo99.res");
  ofstream fConvTSITimeTo99("../OutputFiles/convTSITimeTo99.res");

  fSobolTumDens    << K << " " << 0 << endl;
  fSobolIntTumDens << K << " " << 0 << endl;
  fSobolTimeTo95   << K << " " << 0 << endl;
  fSobolTimeTo99   << K << " " << 0 << endl;

  for(int k(0); k < K; k++){
    fSobolTumDens    << SI[0][k] << " " << TSI[0][k] << endl;
    fSobolIntTumDens << SI[1][k] << " " << TSI[1][k] << endl;
    fSobolTimeTo95   << SI[2][k] << " " << TSI[2][k] << endl;
    fSobolTimeTo99   << SI[3][k] << " " << TSI[3][k] << endl;
     
    for(int l(0); l < NConv; l++){
      fConvSITumDens     << SIConv[l][0][k]  << " ";
      fConvTSITumDens    << TSIConv[l][0][k] << " ";
      fConvSIIntTumDens  << SIConv[l][1][k]  << " ";
      fConvTSIIntTumDens << TSIConv[l][1][k] << " ";
      fConvSITimeTo95    << SIConv[l][2][k]  << " ";
      fConvTSITimeTo95   << TSIConv[l][2][k] << " ";
      fConvSITimeTo99    << SIConv[l][3][k]  << " ";
      fConvTSITimeTo99   << TSIConv[l][3][k] << " ";
    }
    
    fConvSITumDens     << endl;
    fConvTSITumDens    << endl;
    fConvSIIntTumDens  << endl;
    fConvTSIIntTumDens << endl;
    fConvSITimeTo95    << endl;
    fConvTSITimeTo95   << endl;
    fConvSITimeTo99    << endl;
    fConvTSITimeTo99   << endl;	
  }
  
  fSobolTumDens.close();
  fSobolIntTumDens.close();
  fSobolTimeTo95.close();
  fSobolTimeTo99.close();	
  fConvSITumDens.close();
  fConvTSITumDens.close();
  fConvSIIntTumDens.close();
  fConvTSIIntTumDens.close();
  fConvSITimeTo95.close();
  fConvTSITimeTo95.close();
  fConvSITimeTo99.close();
  fConvTSITimeTo99.close();

  free2D(SI, nOut);
  free2D(TSI, nOut);
  free3D(SIConv, NConv, nOut);
  free3D(SIConv, NConv, nOut);
}


void sobolToy(const int N, const string nFRefParInt){
  int const K(5), NConv(log(N) / log(2.0)), nOut(1);
  double h[K], x0[K];
  ifstream fRefParInt(nFRefParInt.c_str());
  
  for(int k(0); k < K; k++){
    fRefParInt >> x0[k];
    fRefParInt >> h[k];
    h[k] -= x0[k];
  }
  
  fRefParInt.close();

  double **SI, **TSI;
  double ***SIConv, ***TSIConv;
  
  SI  = alloc2D(nOut, K);
  TSI = alloc2D(nOut, K);
  SIConv  = alloc3D(NConv, nOut, K);
  TSIConv = alloc3D(NConv, nOut, K);

  sobol(K, N, nOut, x0, h,  SI, TSI, SIConv, TSIConv);

  ofstream fSobolY("../OutputFiles/sobolY.res");
  ofstream fConvSIY("../OutputFiles/convSIY.res");
  ofstream fConvTSIY("../OutputFiles/convTSIY.res");

  fSobolY << K << " " << 0 << endl;
  
  for(int k(0); k < K; k++){
    fSobolY << SI[0][k] << " " << TSI[0][k] << endl;

    for(int l(0); l < NConv; l++){
      fConvSIY  << SIConv[l][0][k]  << " ";
      fConvTSIY << TSIConv[l][0][k] << " ";
    }
    
    fConvSIY  << endl;
    fConvTSIY << endl;
  }
  
  fSobolY.close();	
  fConvSIY.close();
  fConvTSIY.close();

  free2D(SI, nOut);
  free2D(TSI, nOut);
  free3D(SIConv, NConv, nOut);
  free3D(TSIConv, NConv, nOut);
}


void toyModel(double *x, double *y){
  //Legendre polynomial of degree d
  /*int d(x[0]);

    switch (d){
    case 1:
    return x[1];

    case 2:
    return 0.5 * (3.0 * x[1] * x[1] - 1.0);

    case 3:
    return 0.5 * (5.0 * pow(x[1], 3) -
    3.0 * x[1]);

    case 4:
    return 0.125 * (35.0 * pow(x[1], 4) -
    30.0 * x[1] * x[1] + 3.0);

    case 5:
    return 0.125 * (63.0 * pow(x[1], 5) -
    70.0 * pow(x[1], 3) +
    15.0 * x[1]);
    }
    return x[0] * x[1];*/

  y[0] = x[0] + 2.0 * x[1] + x[2] * x[2] + x[3] * x[4];
}
