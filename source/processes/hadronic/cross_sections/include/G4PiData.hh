#ifndef G4PiData_h
#define G4PiData_h

// by J.P Wellisch, Sun Sep 15 2002.

#include "g4std/vector"
#include "g4std/algorithm"
#include "globals.hh"

class G4PiData : public G4std::vector<G4std::pair<G4double, G4std::pair<G4double, G4double > > > 
{
  public:
    G4PiData(const G4double * aTotal, const G4double * aInelastic, const G4double * anEnergy, G4int nPoints);
    struct Delete{void operator()(G4PiData * aP){delete aP;} };
    G4bool AppliesTo(G4double kineticEnergy);
    G4double ReactionXSection(G4double kineticEnergy);
    G4double ElasticXSection(G4double kineticEnergy);
    
  private:
    
};

#endif
