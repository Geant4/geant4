#ifndef G4INTRA_NUCLEI_CASCADER_HH
#define G4INTRA_NUCLEI_CASCADER_HH

#include "G4Collider.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4CascadSpecialFunctions.hh"
#include "G4InuclElementaryParticle.hh"

using namespace G4InuclSpecialFunctions;
using namespace G4CascadSpecialFunctions;

class G4IntraNucleiCascader : public G4Collider {

public:

  G4IntraNucleiCascader() {};

  void setElementaryParticleCollider(G4ElementaryParticleCollider* ecollider) {
    theElementaryParticleCollider = ecollider;   
  };
  
  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				  G4InuclParticle* target);

  void setInteractionCase(G4int intcase) { 
    inter_case = intcase; 
  };

private: 

  G4ElementaryParticleCollider* theElementaryParticleCollider;

  G4int inter_case;

  G4bool goodCase(G4double a, G4double z, G4double eexs, G4double ein) const; 

};        

#endif // G4INTRA_NUCLEI_CASCADER_HH 
