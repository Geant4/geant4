#ifndef G4FISSION_CONFIGURATION_HH
#define G4FISSION_CONFIGURATION_HH

#include "globals.hh"
#include "g4std/iostream"

class G4FissionConfiguration {

public:

  G4FissionConfiguration() {
  };

  G4FissionConfiguration(G4double a, 
			 G4double z, 
			 G4double ez, 
			 G4double ek, 
			 G4double ep) 
    : afirst(a), 
    zfirst(z), 
    ezet(ez), 
    ekin(ek), 
    epot(ep) {
};

  void print() {
    G4cout << " new configuration " << G4endl
	   << " a1 " << afirst << " z1 " << zfirst << " ez " << ezet <<
      " ekin " << ekin << " epot " << epot << G4endl;
  };

  G4double afirst;
  G4double zfirst;
  G4double ezet;
  G4double ekin;
  G4double epot;

};        

#endif // G4FISSION_CONFIGURATION_HH 

