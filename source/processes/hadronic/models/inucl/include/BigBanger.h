#ifndef BIG_BANGER_H
#define BIG_BANGER_H

#include "Collider.h"
#include "InuclElementaryParticle.h"
#include "InuclSpecialFunctions.h"

using namespace InuclSpecialFunctions;

class BigBanger : public Collider {

public:

BigBanger() {};

virtual CollisionOutput collide(InuclParticle* bullet,
                                      InuclParticle* target);

private: 

vector<InuclElementaryParticle>  	    
  generateBangInSCM(double etot, double a, double z, double mp,
	                                              double mn) const;

vector<double> generateMomentumModules(double etot, double a, double z,
       double mp, double mn) const; 

double xProbability(double x, int ia) const; 

double maxProbability(double a) const;

double generateX(int ia, double a, double promax) const; 

};        

#endif // BIG_BANGER_H 
