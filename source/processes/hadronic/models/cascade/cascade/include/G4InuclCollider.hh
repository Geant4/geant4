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
#include "G4IntraNucleiCascader.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InteractionCase.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4Analyser.hh"

using namespace G4InuclSpecialFunctions;

class G4InuclCollider : public G4Collider {

public:

  G4InuclCollider();

  G4InuclCollider(G4ElementaryParticleCollider* ecollider,
		  G4IntraNucleiCascader* incascader, 
		  G4NonEquilibriumEvaporator* noeqevaporator,
		  G4EquilibriumEvaporator* eqevaporator,
		  G4Fissioner* fissioner, 
		  G4BigBanger* bigbanger) {

    setElementaryParticleCollider(ecollider);
    setIntraNucleiCascader(incascader,ecollider);
    setNonEquilibriumEvaporator(noeqevaporator);
    setEquilibriumEvaporator(eqevaporator, fissioner, bigbanger);
    setBigBanger(bigbanger);

  };

  void setElementaryParticleCollider(G4ElementaryParticleCollider* ecollider) {

    theElementaryParticleCollider = ecollider;   
  };

  void setIntraNucleiCascader(G4IntraNucleiCascader* incascader,
			      G4ElementaryParticleCollider* ecollider) {

    theIntraNucleiCascader = incascader;
    theIntraNucleiCascader->setElementaryParticleCollider(ecollider);
  };

  void setNonEquilibriumEvaporator(G4NonEquilibriumEvaporator* noeqevaporator) {

    theNonEquilibriumEvaporator = noeqevaporator;   
  };

  void setEquilibriumEvaporator(G4EquilibriumEvaporator* eqevaporator,
				G4Fissioner* fissioner, 
				G4BigBanger* bigbanger) {

    theEquilibriumEvaporator = eqevaporator;
    theEquilibriumEvaporator->setFissioner(fissioner);   
    theEquilibriumEvaporator->setBigBanger(bigbanger);   
  };

  void setBigBanger(G4BigBanger* bigbanger) {

    theBigBanger = bigbanger;   
  };
  
  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target);

private: 

  G4int verboseLevel;

  G4bool inelasticInteractionPossible(G4InuclParticle* bullet,
				      G4InuclParticle* target, 
				      G4double ekin) const;

  G4InteractionCase bulletTargetSetter(G4InuclParticle* bullet,
				       G4InuclParticle* target) const; 

  G4bool explosion(G4InuclNuclei* target) const; 
       
  G4ElementaryParticleCollider* theElementaryParticleCollider;
  G4IntraNucleiCascader* theIntraNucleiCascader;
  G4NonEquilibriumEvaporator* theNonEquilibriumEvaporator;
  G4EquilibriumEvaporator* theEquilibriumEvaporator;
  G4BigBanger* theBigBanger;

};        

#endif // G4INUCL_COLLIDER_HH 


