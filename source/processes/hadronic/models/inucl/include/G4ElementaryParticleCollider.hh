#ifndef G4ELEMENTARY_PARTICLE_COLLIDER_HH
#define G4ELEMENTARY_PARTICLE_COLLIDER_HH

#include "G4Collider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4CascadSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"

using namespace G4InuclSpecialFunctions;
using namespace G4CascadSpecialFunctions;

class G4ElementaryParticleCollider : public G4Collider {

public:

  G4ElementaryParticleCollider() {};

  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target);

private: 

  G4int generateMultiplicity(G4int is, 
			     G4double ekin) const;
      
  vector<G4InuclElementaryParticle> generateSCMfinalState(G4double ekin, 
							  G4double etot_scm, G4double pscm,	     
							  G4InuclElementaryParticle* particle1,
							  G4InuclElementaryParticle* particle2, 
							  G4LorentzConvertor* toSCM) const; 

  vector<G4double> generateMomModules(const vector<G4int>& kinds, 
				      G4int mult,
				      G4int is, 
				      G4double ekin, 
				      G4double etot_cm) const; 
      
  G4bool reChargering(G4double ekin, 
		      G4int is) const;

  vector<G4double> particleSCMmomentumFor2to2(G4int is, 
					      G4int kw, 
					      G4double ekin,
					      G4double pscm) const; 
	    
  G4int getElasticCase(G4int is, 
		       G4int kw, 
		       G4double ekin) const;

  vector<G4int> generateOutgoingKindsFor2toMany(G4int is, 
						G4int mult, 
						G4double ekin) const; 

  G4double getMomModuleFor2toMany(G4int is, 
				  G4int mult, 
				  G4int knd, 
				  G4double ekin) const; 

  G4bool satisfyTriangle(const vector<G4double>& modules) const; 
	
  vector<G4double> particleSCMmomentumFor2to3(G4int is, 
					      G4int knd, 
					      G4double ekin, 
					      G4double pmod) const; 
	
  G4int getIL(G4int is, 
	      G4int mult) const; 

  pair<G4double, G4double> adjustIntervalForElastic(G4double ekin, 
						    G4double ak, 
						    G4double ae,
						    G4int k, 
						    G4int l, 
						    const vector<G4double>& ssv, 
						    G4double st) const;
 
  vector<G4InuclElementaryParticle> 
  generateSCMpionAbsorption(G4double etot_scm,
			    G4InuclElementaryParticle* particle1,
			    G4InuclElementaryParticle* particle2) const; 
    
};        

#endif // G4ELEMENTARY_PARTICLE_COLLIDER_HH 
















