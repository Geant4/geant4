#ifndef FISSIONER_H
#define FISSIONER_H

#include "Collider.h"
#include "InuclSpecialFunctions.h"

using namespace InuclSpecialFunctions;

class Fissioner : public Collider {

public:

Fissioner() {};

virtual CollisionOutput collide(InuclParticle* bullet,
                     InuclParticle* target);

private: 

double getC2(double A1, double A2, double X3, double X4, double R12) const; 

double getZopt(double A1, double A2, double ZT, 
                   double X3, double X4, double R12) const;
		    
void potentialMinimization(double& VP, vector<double>& ED, double& VC,
   double AF, double AS, double ZF, double ZS,
   vector<double>& AL1, vector<double>& BET1, double& R12) const; 

};        

#endif // FISSIONER_H 
