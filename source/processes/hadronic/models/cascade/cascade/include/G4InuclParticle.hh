//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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








