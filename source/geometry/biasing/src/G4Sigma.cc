#include "G4Sigma.hh"
#include "globals.hh"
#include <cstdio>
#include <cmath>
 
G4Sigma::G4Sigma(){
  Init();
}

G4Sigma::~G4Sigma(){}

void G4Sigma::Init() {
  fEntries = 0;
  fMean = 0.;
  fSigma = -1.;
  fXsum = 0;
  fXXsum = 0;
  fWsum = 0;
  fWXsum = 0;
  fWXXsum = 0;
  fcalc = 0;
}

int G4Sigma::Xin(double x, double w){
  if (w<0.) Error("Xin: w < 0");
  fEntries++;
  fXsum+=x;
  fXXsum+=x*x;
  fWsum += w;
  fWXsum += w*x;
  fWXXsum += w*x*x;
  fcalc = 0;
  return fEntries;
}



int G4Sigma::GetCalc() const {
  return fcalc;
}

int G4Sigma::Calculate() const {
  if (fcalc==0) {
    if(fWsum>0) {
      fMean=fWXsum/fWsum;
      fSigma = sqrt( fWXXsum / fWsum - fMean * fMean);
      fcalc = 1;
    } else {
      fcalc = -1;
    }
  }
  return fcalc;
}


int G4Sigma::GetEntries() const {
  return fEntries;
}

double G4Sigma::GetMean() const {
  if (fcalc==0) Calculate();
  return fMean;
}

double G4Sigma::GetSigma() const {
  if (fcalc==0) Calculate();
  return fSigma;
}

double G4Sigma::GetXsum() const {return fXsum;}
double G4Sigma::GetXXsum() const {return fXXsum;}
double G4Sigma::GetSumOfWeights() const {return fWsum;}
double G4Sigma::GetWeightedXsum() const {return fWXsum;}
double G4Sigma::GetWeightedXXsum() const {return fWXXsum;}
void G4Sigma::Error(G4String m){
  G4cout << "Error:G4Sigma::" << m << G4endl;
}


ostream& operator<<(ostream &out, const G4Sigma &s){
  out << "entries                             : " << s.GetEntries() << "\n";
  out << "Sum(w)                              : " << s.GetSumOfWeights()<<"\n";
  out << "Sum(w*x)                            : " << s.GetWeightedXsum() << "\n";
  out << "Sum(w*x*x)                          : " << s.GetWeightedXXsum() << "\n";
  out << "mean=Sum(w*x) / Sum(w)              : " << s.GetMean() << "\n";
  out << "sigma=sqrt(Sum(w*x*x)/Sum(w)-mean^2): " << s.GetSigma() << "\n";
  out << "Sum(x)                              : " << s.GetXsum() << "\n";
  out << "Sum(x^2)                            : " << s.GetXXsum() << "\n";
  
  return out;
};



