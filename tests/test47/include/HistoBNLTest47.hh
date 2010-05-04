#ifndef HistoBNLTest47_H
#define HistoBNLTest47_H

#include "HistoTest47.hh"
#include "HistoEPTest47.hh"

#include "TFile.h"
#include "TH1F.h"
#include <string>

class HistoBNLTest47 : public HistoTest47 {

public:

  HistoBNLTest47(std::string namePart, std::string nameMat, G4double momentum,
		 std::string nameGen);
  virtual ~HistoBNLTest47();

  virtual void fill(G4VParticleChange*, G4LorentzVector);
  virtual void write(G4double cross_sec, G4int nevt) ;

private:

  void initialize();
  void book();

  std::string           fileName;
  char                  tag1Name[60], tag2Name[24], tag3Name[40];
  HistoEPTest47         epTest;
  std::vector<TH1F*>    hiMT11, hiMT12, hiMT21, hiMT22;
  std::vector<TH1F*>    hiMT31, hiMT32, hiMT41, hiMT42;
  std::vector<TH1F*>    hiMT51, hiMT52, hiMT61, hiMT62;
  std::vector<G4double> rapidities, ymin, ymax;
  
};

#endif
