#ifndef ANAParticle_h
#define ANAParticle_h

#include "globals.hh"
#include <fstream>

class ANAParticle
{
  public:
  
    G4bool Init(ifstream & theData);
    
    G4int GetCharge(){return charge;}
    G4int GetPDGCode(){return pdgCode;}
    G4double GetEnergy() {return energy;}
    G4double GetKineticEnergy()
    {
      return energy - GetMomentum();
    }
    G4double GetPx() {return px;}
    G4double GetPy() {return py;}
    G4double GetPz() {return pz;}
    G4double GetMomentum()
    {
      G4double result = sqrt(max(G4double(0), px*px+py*py+pz*pz));
      return result;
    }
    G4double GetCosTheta() { return pz/GetMomentum(); }
    
    G4double GetWeight() 
    {
      G4double result = 1./(4.*3.14159265*GetMomentum());
      return result;
    }
    
  private:
  G4int charge;
  G4int pdgCode;
  G4double energy;
  G4double px;
  G4double py;
  G4double pz;
};

inline
G4bool ANAParticle::Init(ifstream & theData)
{
  G4int dummy;
  if(!(theData>>charge)) return false;
  theData>>pdgCode;
  theData>>energy;
  theData>>px;
  theData>>py;
  theData>>pz;
  theData>>dummy;
  return true;
}

#endif
