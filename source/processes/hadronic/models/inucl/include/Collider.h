#ifndef COLLIDER_H
#define COLLIDER_H

#include "InuclParticle.h"
#include "CollisionOutput.h"

class Collider {

public:

Collider() {};

virtual CollisionOutput collide(InuclParticle* bullet,
                                       InuclParticle* target) = 0;

};        

#endif // COLLIDER_H 
