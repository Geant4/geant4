
#ifndef __NewElastic__
#define __NewElastic__


#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4Track.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4LightMedia.hh"


class HadronElastic : public G4HadronicInteraction
{
public:
  HadronElastic() : G4HadronicInteraction()
  {
    SetMinEnergy(0*MeV);
    SetMaxEnergy(10*GeV);
    CalcComb();
    //    PrintComb();
  }
  ~HadronElastic(){};
  G4VParticleChange* ApplyYourself(const G4Track& aTrack, G4Nucleus& targetNucleus);
private:
  G4double GetMomentumHint(G4double num,G4double maxP,G4Nucleus& targetNucl,G4int Z);
  G4double GetMaximalValue(G4double maxP,G4Nucleus& targetNucl);
  G4LightMedia LightMedia;
  void CalcComb(void);
  void PrintComb();
  G4double GetComb(unsigned A,unsigned n);
  G4double m_Comb[20100];
};

#endif
