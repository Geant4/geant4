#ifndef G4Sigma_hh
#define G4Sigma_hh 1
#include "globals.hh"
#include <iostream>


class G4Sigma { 
  
 public:
  G4Sigma();
  ~G4Sigma();
  void Init();
  int Xin(double x, double w = 1.);
  int GetEntries() const;
  double GetMean() const;
  double GetSigma() const;
  double GetXsum() const;
  double GetXXsum() const;
  double GetSumOfWeights() const;
  double GetWeightedXsum() const;
  double GetWeightedXXsum() const;

 private:
  int GetCalc() const;
  int Calculate() const;
  int fEntries;
  mutable double fMean;
  mutable double fSigma;
  double fXsum;
  double fXXsum;
  double fWsum;
  double fWXsum;
  double fWXXsum;
  mutable int fcalc;
  void Error(G4String m);

};

G4std::ostream& operator<<(G4std::ostream &out, const G4Sigma &s);

#endif




