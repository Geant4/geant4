#ifndef G4INUCL_NUCLEI_HH
#define G4INUCL_NUCLEI_HH

#include "G4InuclParticle.hh"
#include "G4ExitonConfiguration.hh"
#include "G4InuclSpecialFunctions.hh"

using namespace G4InuclSpecialFunctions;

class G4InuclNuclei : public G4InuclParticle {

public:

G4InuclNuclei() {};

G4InuclNuclei(G4double a, G4double z) : A(a), Z(z) { setNucleiMass();
     exitationEnergy = 0.0; };

G4InuclNuclei(const vector<G4double>& mom, G4double a, G4double z) :
        InuclParticle(mom), A(a), Z(z) { setNucleiMass();  
	exitationEnergy = 0.0;};

G4InuclNuclei(G4double ekin, G4double a, G4double z) : A(a), Z(z) { 
  setNucleiMass();
  vector<G4double> mom(4, 0.0);
  mom[0] = ekin + nucleiMass;
  mom[3] = sqrt(mom[0]*mom[0] - nucleiMass*nucleiMass);
  G4InuclParticle::setMomentum(mom);
  exitationEnergy = 0.0;
};
	
void setA(G4double a) { A = a; };

void setZ(G4double z) { Z = z; };

void setExitationEnergy(G4double e) { exitationEnergy = e; };

void setEnergy() {
  momentum[0] = sqrt(momentum[1]*momentum[1] + momentum[2]*momentum[2] +
     momentum[3]*momentum[3] + nucleiMass*nucleiMass);
  
};

void setNucleiMass() {
  nucleiMass = 0.93827*Z + 0.93957*(A - Z) - 0.001*bindingEnergy(A,Z);
};

void setExitonConfiguration(const G4ExitonConfiguration& config) { 
  theExitonConfiguration = config;
};

G4double getA() const { return A; };

G4double getZ() const { return Z; };

G4double getExitationEnergy() const { return exitationEnergy; };

G4double getExitationEnergyInGeV() const { return 0.001*exitationEnergy; };

G4ExitonConfiguration getExitonConfiguration() const {
  return theExitonConfiguration;
};

G4double getMass() const { return nucleiMass; };

G4double getKineticEnergy() const { return momentum[0] - nucleiMass; };

G4double getNucleiMass(G4double a, G4double z) const {
  return 0.93827*z + 0.93957*(a - z) - 0.001*bindingEnergy(a, z);
};

virtual void printParticle() const {
  G4cout << " A " << A << " Z " << Z << " mass " << nucleiMass << 
       " Eex (MeV) " << exitationEnergy << G4endl;
  if(momentum.size() == 4)         
    G4cout << " Px " << momentum[1] << " Py " << momentum[2] << " Pz " <<  
          momentum[3] <<  " E " << momentum[0] << G4endl;
};

private: 

G4double A;
G4double Z;
G4double exitationEnergy;
G4double nucleiMass;
G4ExitonConfiguration theExitonConfiguration;

};        

#endif // G4INUCL_NUCLEI_HH 






