#ifndef G4NON_EQUILIBRIUM_EVAPORATOR_HH
#define G4NON_EQUILIBRIUM_EVAPORATOR_HH

#include "G4Collider.hh"
#include "G4InuclSpecialFunctions.hh"

using namespace G4InuclSpecialFunctions;

class G4NonEquilibriumEvaporator : public G4Collider {

public:

  G4NonEquilibriumEvaporator();

  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target);

private: 
G4int verboseLevel;
  G4double getMatrixElement(G4double A) const;

  G4double getE0(G4double A) const; 

  G4double getParLev(G4double A, G4double Z) const;

};

#endif // G4NON_EQUILIBRIUM_EVAPORATOR_HH 
