#ifndef EQUILIBRIUM_EVAPORATOR_H
#define EQUILIBRIUM_EVAPORATOR_H

#include "Collider.h"
#include "Fissioner.h"
#include "BigBanger.h"
#include "InuclSpecialFunctions.h"

using namespace InuclSpecialFunctions;

class EquilibriumEvaporator : public Collider {

public:

EquilibriumEvaporator() {};

void setFissioner(Fissioner* fissioner) {
    theFissioner = fissioner;
};

void setBigBanger(BigBanger* banger) {
    theBigBanger = banger;
};

virtual CollisionOutput collide(InuclParticle* bullet,
                                   InuclParticle* target);

private: 

double getE0(double A) const; 

double getPARLEVDEN(double A, double Z) const; 

bool timeToBigBang(double a, double z, double e) const;

bool goodRemnant(double a, double z) const; 

double getQF(double x, double x2, double a, double z, double e) const;

double getAF(double x, double a, double z, double e) const; 

Fissioner* theFissioner;
BigBanger* theBigBanger;

};        

#endif // EQUILIBRIUM_EVAPORATOR_H 
