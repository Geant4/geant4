#ifndef ELEMENTARY_PARTICLE_COLLIDER_H
#define ELEMENTARY_PARTICLE_COLLIDER_H

#include "Collider.h"
#include "InuclElementaryParticle.h"
#include "InuclSpecialFunctions.h"
#include "CascadSpecialFunctions.h"
#include "LorentzConvertor.h"

using namespace InuclSpecialFunctions;
using namespace CascadSpecialFunctions;

class ElementaryParticleCollider : public Collider {

public:

ElementaryParticleCollider() {};

virtual CollisionOutput collide(InuclParticle* bullet,
                     InuclParticle* target);

private: 

int generateMultiplicity(int is, double ekin) const;
      
vector<InuclElementaryParticle> generateSCMfinalState(double ekin, 
       double etot_scm, double pscm,	     
      InuclElementaryParticle* particle1,
      InuclElementaryParticle* particle2, LorentzConvertor* toSCM) const; 

vector<double> 	generateMomModules(const vector<int>& kinds, int mult,
      int is, double ekin, double etot_cm) const; 
      
bool reChargering(double ekin, int is) const;

vector<double> particleSCMmomentumFor2to2(int is, int kw, double ekin,
            double pscm) const; 
	    
int getElasticCase(int is, int kw, double ekin) const;

vector<int> generateOutgoingKindsFor2toMany(
        int is, int mult, double ekin) const; 

double getMomModuleFor2toMany( 
        int is, int mult, int knd, double ekin) const; 

bool satisfyTriangle(const vector<double>& modules) const; 
	
vector<double> particleSCMmomentumFor2to3(
        int is, int knd, double ekin, double pmod) const; 
	
int getIL(int is, int mult) const; 

pair<double,double> adjustIntervalForElastic(double ekin, double ak, double ae,
 int k, int l, const vector<double>& ssv, double st) const;
 
vector<InuclElementaryParticle> 
    generateSCMpionAbsorption(double etot_scm,
      InuclElementaryParticle* particle1,
      InuclElementaryParticle* particle2) const; 
    
};        

#endif // ELEMENTARY_PARTICLE_COLLIDER_H 
