#ifndef G4INUCL_PARTICLE_HH
#define G4INUCL_PARTICLE_HH

#ifndef GLOB
#include "globals.hh"
#endif

#include "g4std/iostream"
#include "g4std/vector"

class G4InuclParticle {

public:

  G4InuclParticle() {};

  virtual ~G4InuclParticle() { };
 
  G4InuclParticle(const G4std::vector<G4double>& mom) {

    setMomentum(mom);
  };

  void setMomentum(const G4std::vector<G4double>& mom) {

    momentum = mom;
  };

  G4std::vector<G4double> getMomentum() const { 

    return momentum; 
  };

  G4double getMomModule() const { 

    return sqrt(momentum[1] * momentum[1] +
		momentum[2] * momentum[2] + 
		momentum[3] * momentum[3]); 
  };
   
  virtual void printParticle() const {

    G4cout << " px " << momentum[1] << " py " << momentum[2] <<
      " pz " << momentum[3] <<
      " pmod " << sqrt(momentum[1] * momentum[1] + 
		       momentum[2] * momentum[2] +
		       momentum[3] * momentum[3])
	   << " E " << momentum[0] << G4endl;
  };

protected: 

  G4std::vector<G4double> momentum;

};        

#endif // G4INUCL_PARTICLE_HH 








