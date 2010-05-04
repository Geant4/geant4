#ifndef HistoEPTest47_H
#define HistoEPTest47_H

#include "HistoTest47.hh"

#include "TFile.h"
#include "TH1F.h"

class HistoEPTest47 {

public:

  HistoEPTest47();
  virtual ~HistoEPTest47();

  void initialize(std::string namePart, std::string nameMat, G4double momentum,
		  std::string nameGen);
  void fill(G4VParticleChange*, G4LorentzVector);
  void write() ;

private:

  TH1F                  *hiPx, *hiPy, *hiPz, *hiE, *hiKE;
 
};

#endif
