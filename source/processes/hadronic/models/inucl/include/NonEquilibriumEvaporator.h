#ifndef NON_EQUILIBRIUM_EVAPORATOR_H
#define NON_EQUILIBRIUM_EVAPORATOR_H

#include "Collider.h"
#include "InuclSpecialFunctions.h"

using namespace InuclSpecialFunctions;

class NonEquilibriumEvaporator : public Collider {

public:

NonEquilibriumEvaporator() {};

virtual CollisionOutput collide(InuclParticle* bullet,
                     InuclParticle* target);

private: 

double getMatrixElement(double A) const;

double getE0(double A) const; 

double getParLev(double A, double Z) const;

};

#endif // NON_EQUILIBRIUM_EVAPORATOR_H 
