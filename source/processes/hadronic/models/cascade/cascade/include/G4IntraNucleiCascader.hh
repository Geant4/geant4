//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#include "G4Collider.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4CascadSpecialFunctions.hh"
#include "G4InuclElementaryParticle.hh"

using namespace G4InuclSpecialFunctions;
using namespace G4CascadSpecialFunctions;

class G4IntraNucleiCascader : public G4Collider {

public:

  G4IntraNucleiCascader();

  void setElementaryParticleCollider(G4ElementaryParticleCollider* ecollider) {
    theElementaryParticleCollider = ecollider;   
  };
  
  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				  G4InuclParticle* target);

  void setInteractionCase(G4int intcase) { 
    inter_case = intcase; 
  };

private: 
G4int verboseLevel;
  G4ElementaryParticleCollider* theElementaryParticleCollider;

  G4int inter_case;

  G4bool goodCase(G4double a, G4double z, G4double eexs, G4double ein) const; 

};        

#endif // G4INTRA_NUCLEI_CASCADER_HH 
