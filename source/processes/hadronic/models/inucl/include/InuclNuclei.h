#ifndef INUCL_NUCLEI_H
#define INUCL_NUCLEI_H

#include "InuclParticle.h"
#include "ExitonConfiguration.h"
#include "InuclSpecialFunctions.h"

using namespace InuclSpecialFunctions;

class InuclNuclei : public InuclParticle {

public:

InuclNuclei() {};

InuclNuclei(double a, double z) : A(a), Z(z) { setNucleiMass();
     exitationEnergy = 0.; };

InuclNuclei(const vector<double>& mom, double a, double z) :
        InuclParticle(mom), A(a), Z(z) { setNucleiMass();  
	exitationEnergy = 0.;};

InuclNuclei(double ekin, double a, double z) : A(a), Z(z) { 
  setNucleiMass();
  vector<double> mom(4,0.);
  mom[0] = ekin + nucleiMass;
  mom[3] = sqrt(mom[0]*mom[0] - nucleiMass*nucleiMass);
  InuclParticle::setMomentum(mom);
  exitationEnergy = 0.;
};
	
void setA(double a) { A = a; };

void setZ(double z) { Z = z; };

void setExitationEnergy(double e) { exitationEnergy = e; };

void setEnergy() {
  momentum[0] = sqrt(momentum[1]*momentum[1] + momentum[2]*momentum[2] +
     momentum[3]*momentum[3] + nucleiMass*nucleiMass);
  
};

void setNucleiMass() {
  nucleiMass = 0.93827*Z + 0.93957*(A - Z) - 0.001*bindingEnergy(A,Z);
};

void setExitonConfiguration(const ExitonConfiguration& config) { 
  theExitonConfiguration = config;
};

double getA() const { return A; };

double getZ() const { return Z; };

double getExitationEnergy() const { return exitationEnergy; };

double getExitationEnergyInGeV() const { return 0.001*exitationEnergy; };

ExitonConfiguration getExitonConfiguration() const {
  return theExitonConfiguration;
};

double getMass() const { return nucleiMass; };

double getKineticEnergy() const { return momentum[0] - nucleiMass; };

double getNucleiMass(double a, double z) const {
  return 0.93827*z + 0.93957*(a - z) - 0.001*bindingEnergy(a,z);
};

virtual void printParticle() const {
  cout << " A " << A << " Z " << Z << " mass " << nucleiMass << 
       " Eex (MeV) " << exitationEnergy << endl;
  if(momentum.size() == 4)         
    cout << " Px " << momentum[1] << " Py " << momentum[2] << " Pz " <<  
          momentum[3] <<  " E " << momentum[0] << endl;
};

private: 

double A;
double Z;
double exitationEnergy;
double nucleiMass;
ExitonConfiguration theExitonConfiguration;

};        

#endif // INUCL_NUCLEI_H 
