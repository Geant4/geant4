#ifndef INUCL_COLLIDER_H
#define INUCL_COLLIDER_H

#include "Collider.h"
#include "IntraNucleiCascader.h"
#include "NonEquilibriumEvaporator.h"
#include "EquilibriumEvaporator.h"
#include "Fissioner.h"
#include "BigBanger.h"
#include "ElementaryParticleCollider.h"
#include "InteractionCase.h"
#include "InuclNuclei.h"
#include "InuclSpecialFunctions.h"
#include "Analyser.h"

using namespace InuclSpecialFunctions;

class InuclCollider : public Collider {

public:

InuclCollider() {};

InuclCollider(ElementaryParticleCollider* ecollider,
              IntraNucleiCascader* incascader, 
              NonEquilibriumEvaporator* noeqevaporator,
	      EquilibriumEvaporator* eqevaporator,
	      Fissioner* fissioner, BigBanger* bigbanger) {
  setElementaryParticleCollider(ecollider);
  setIntraNucleiCascader(incascader,ecollider);
  setNonEquilibriumEvaporator(noeqevaporator);
  setEquilibriumEvaporator(eqevaporator,fissioner,bigbanger);
  setBigBanger(bigbanger);
};

void setElementaryParticleCollider(ElementaryParticleCollider* ecollider) {
  theElementaryParticleCollider = ecollider;   
};

void setIntraNucleiCascader(IntraNucleiCascader* incascader,
  ElementaryParticleCollider* ecollider) {
  theIntraNucleiCascader = incascader;
  theIntraNucleiCascader->setElementaryParticleCollider(ecollider);
};

void setNonEquilibriumEvaporator(NonEquilibriumEvaporator* noeqevaporator) {
  theNonEquilibriumEvaporator = noeqevaporator;   
};

void setEquilibriumEvaporator(EquilibriumEvaporator* eqevaporator,
    Fissioner* fissioner, BigBanger* bigbanger) {
  theEquilibriumEvaporator = eqevaporator;
  theEquilibriumEvaporator->setFissioner(fissioner);   
  theEquilibriumEvaporator->setBigBanger(bigbanger);   
};

void setBigBanger(BigBanger* bigbanger) {
  theBigBanger = bigbanger;   
};
  
virtual CollisionOutput collide(InuclParticle* bullet,
                     InuclParticle* target);

private: 

bool inelasticInteractionPossible(InuclParticle* bullet,
          InuclParticle* target, double ekin) const;

InteractionCase bulletTargetSetter(InuclParticle* bullet,
          InuclParticle* target) const; 

bool explosion(InuclNuclei* target) const; 
       
ElementaryParticleCollider* theElementaryParticleCollider;
IntraNucleiCascader* theIntraNucleiCascader;
NonEquilibriumEvaporator* theNonEquilibriumEvaporator;
EquilibriumEvaporator* theEquilibriumEvaporator;
BigBanger* theBigBanger;

};        

#endif // INUCL_COLLIDER_H 
