#ifndef NUCLEI_MODEL_H
#define NUCLEI_MODEL_H

#include "InuclNuclei.h"
#include "InuclElementaryParticle.h"
#include "CascadParticle.h"
#include "InuclSpecialFunctions.h"
#include "CascadSpecialFunctions.h"
#include "InuclElementaryParticle.h"
#include "ElementaryParticleCollider.h"
#include "vector"

using namespace InuclSpecialFunctions;
using namespace CascadSpecialFunctions;

typedef pair<InuclElementaryParticle,double> partner;
typedef vector<partner> partners;

class NucleiModel {

public:

NucleiModel() {};

NucleiModel(InuclNuclei* nuclei) {
  generateModel(nuclei->getA(), nuclei->getZ());
};

void generateModel(double a, double z);

void reset() {
  neutronNumberCurrent = neutronNumber;
  protonNumberCurrent = protonNumber;
};

void printModel() const; 

double getDensity(int ip, int izone) const {
  return nucleon_densities[ip-1][izone];
};

double getFermiMomentum(int ip, int izone) const {
  return fermi_momenta[ip-1][izone];
};

double getFermiKinetic(int ip, int izone) const {
  double ekin = 0.;
  if(ip < 3 && izone < number_of_zones) {
    double pf = fermi_momenta[ip-1][izone]; 
    double mass = ip == 1 ? 0.93827 : 0.93957;
    ekin = sqrt(pf*pf + mass*mass) - mass;
  };  
  return ekin;
};

double getPotential(int ip, int izone) const {
  int ip0 = ip < 3 ? ip - 1 : 2;
  return izone < number_of_zones ? zone_potentials[ip0][izone] : 0.;
};

vector<CascadParticle> generateParticleFate(CascadParticle& cparticle,
     ElementaryParticleCollider* theElementaryParticleCollider); 

double getNumberOfNeutrons() const { return neutronNumberCurrent; };

double getNumberOfProtons() const { return protonNumberCurrent; };

bool empty() const { return neutronNumberCurrent < 1. && 
     protonNumberCurrent < 1.; };

bool stillInside(const CascadParticle& cparticle) {
  return cparticle.getCurrentZone() < number_of_zones;
};

CascadParticle initializeCascad(InuclElementaryParticle* particle);

pair<vector<CascadParticle>, vector<InuclElementaryParticle> >
   NucleiModel::initializeCascad(InuclNuclei* bullet, InuclNuclei* target);

pair<int,int> getTypesOfNucleonsInvolved() const {
  return pair<int,int>(current_nucl1,current_nucl2);
};
bool worthToPropagate(const CascadParticle& cparticle) const; 
    
private: 

bool passFermi(const vector<InuclElementaryParticle>& particles, int zone);

void boundaryTransition(CascadParticle& cparticle);

InuclElementaryParticle generateNucleon(int type, int zone) const;

InuclElementaryParticle generateQuasiDeutron(int type1, int type2,
                   int zone) const;

partners generateInteractionPartners(CascadParticle& cparticle) const;

double volNumInt(double r1, double r2, double cu, double d1) const; 

double volNumInt1(double r1, double r2, double cu2) const; 

double getRatio(int ip) const;

vector<vector<double> > nucleon_densities;

vector<vector<double> > zone_potentials;

vector<vector<double> > fermi_momenta;

vector<double> zone_radii;

vector<double> binding_energies;

double nuclei_radius;

int number_of_zones;

double A;
double Z;
double neutronNumber;
double protonNumber;
double neutronNumberCurrent;
double protonNumberCurrent;

int current_nucl1;
int current_nucl2;
 
};        

#endif // NUCLEI_MODEL_H 
