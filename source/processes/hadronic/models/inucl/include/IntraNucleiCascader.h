#ifndef INTRA_NUCLEI_CASCADER_H
#define INTRA_NUCLEI_CASCADER_H

#include "Collider.h"
#include "ElementaryParticleCollider.h"
#include "InuclSpecialFunctions.h"
#include "CascadSpecialFunctions.h"
#include "InuclElementaryParticle.h"

using namespace InuclSpecialFunctions;
using namespace CascadSpecialFunctions;

class IntraNucleiCascader : public Collider {

public:

IntraNucleiCascader() {};

void setElementaryParticleCollider(ElementaryParticleCollider* ecollider) {
  theElementaryParticleCollider = ecollider;   
};
  
virtual CollisionOutput collide(InuclParticle* bullet,
                     InuclParticle* target);

void setInteractionCase(int intcase) { inter_case = intcase; };

private: 

ElementaryParticleCollider* theElementaryParticleCollider;

int inter_case;

bool goodCase(double a, double z, double eexs, double ein) const; 

};        

#endif // INTRA_NUCLEI_CASCADER_H 
