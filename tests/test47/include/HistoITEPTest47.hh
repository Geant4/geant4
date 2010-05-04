#ifndef HistoITEPTest47_H
#define HistoITEPTest47_H

#include "HistoTest47.hh"
#include "HistoEPTest47.hh"

#include "TFile.h"
#include "TH1F.h"
#include <string>

class G4VParticleChange;

class HistoITEPTest47 : public HistoTest47 {

public:

  HistoITEPTest47(std::string namePart, std::string nameMat, G4double momentum,
		  std::string nameGen);
  virtual ~HistoITEPTest47();

  virtual void fill(G4VParticleChange*, G4LorentzVector);
  virtual void write(G4double cross_sec, G4int nevt);

private:

  void initialize();
  void book();

  std::string           fileName;
  char                  tag1Name[60], tag2Name[24], tag3Name[40];
  double                dtheta, de;
  HistoEPTest47         epTest;
  std::vector<TH1F*>    hiKE11, hiKE12, hiCT11, hiCT12;
  std::vector<TH1F*>    hiKE21, hiKE22, hiCT21, hiCT22;
  std::vector<G4double> energies, emin, emax, angles, cthmin, cthmax, dcth;
  
};

#endif
