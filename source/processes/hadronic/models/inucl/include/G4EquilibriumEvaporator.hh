#ifndef G4EQUILIBRIUM_EVAPORATOR_HH
#define G4EQUILIBRIUM_EVAPORATOR_HH

#include "G4Collider.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4InuclSpecialFunctions.hh"

using namespace G4InuclSpecialFunctions;

class G4EquilibriumEvaporator : public G4Collider {

public:

  G4EquilibriumEvaporator() {};

  void setFissioner(G4Fissioner* fissioner) {
    theFissioner = fissioner;
  };

  void setBigBanger(G4BigBanger* banger) {
    theBigBanger = banger;
  };

  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target);

private: 

  G4double getE0(G4double A) const; 

  G4double getPARLEVDEN(G4double A, 
			G4double Z) const; 

  G4bool timeToBigBang(G4double a, 
		       G4double z, 
		       G4double e) const;

  G4bool goodRemnant(G4double a, 
		     G4double z) const; 

  G4double getQF(G4double x, 
		 G4double x2, 
		 G4double a, 
		 G4double z, 
		 G4double e) const;

  G4double getAF(G4double x, 
		 G4double a, 
		 G4double z, 
		 G4double e) const; 

  G4Fissioner* theFissioner;
  G4BigBanger* theBigBanger;

};        

#endif // G4EQUILIBRIUM_EVAPORATOR_HH 


