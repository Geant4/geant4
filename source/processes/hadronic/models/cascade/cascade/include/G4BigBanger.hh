#ifndef G4BIG_BANGER_HH
#define G4BIG_BANGER_HH

#include "G4Collider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclSpecialFunctions.hh"


using namespace G4InuclSpecialFunctions;

class G4BigBanger : public G4Collider {

public:

  G4BigBanger();

  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target);

private: 

G4int verboseLevel;
  vector<G4InuclElementaryParticle> generateBangInSCM(G4double etot, 
						      G4double a, 
						      G4double z, 
						      G4double mp,
						      G4double mn) const;

  vector<G4double> generateMomentumModules(G4double etot, 
					   G4double a, 
					   G4double z,
					   G4double mp, 
					   G4double mn) const; 

  G4double xProbability(G4double x, 
			G4int ia) const; 

  G4double maxProbability(G4double a) const;

  G4double generateX(G4int ia, 
		     G4double a, 
		     G4double promax) const; 

};        

#endif // G4BIG_BANGER_HH 











