//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef G4INUCL_NUCLEI_HH
#define G4INUCL_NUCLEI_HH

#ifndef G4INUCL_PARTICLE_HH
#include "G4InuclParticle.hh"
#endif
#include "G4ExitonConfiguration.hh"
#include "G4InuclSpecialFunctions.hh"

using namespace G4InuclSpecialFunctions;

class G4InuclNuclei : public G4InuclParticle {

public:

  G4InuclNuclei() {};

  G4InuclNuclei(G4double a, 
		G4double z)
    : A(a),
      Z(z) { 

    setNucleiMass();
    exitationEnergy = 0.0; 
  };

  G4InuclNuclei(const G4CascadeMomentum& mom, 
		G4double a, 
		G4double z)
    : G4InuclParticle(mom),
      A(a),
      Z(z) { 

    setNucleiMass();  
    exitationEnergy = 0.0;
  };

  G4InuclNuclei(G4double ekin, 
		G4double a, 
		G4double z) 
    : A(a), 
    Z(z) { 

    setNucleiMass();
    G4CascadeMomentum mom;
    mom[0] = ekin + nucleiMass;
    mom[3] = std::sqrt(mom[0] * mom[0] - nucleiMass * nucleiMass);
    G4InuclParticle::setMomentum(mom);
    exitationEnergy = 0.0;
  };
	
  void setA(G4double a) { 

    A = a; 
  };

  void setZ(G4double z) { 

    Z = z; 
  };

  void setExitationEnergy(G4double e) { 

    exitationEnergy = e; 
  };

  void setEnergy() {

    momentum[0] = std::sqrt(momentum[1] * momentum[1] + momentum[2] * momentum[2] +
		       momentum[3] * momentum[3] + nucleiMass * nucleiMass);  
  };

  void setNucleiMass() {

    nucleiMass = 0.93827 * Z + 0.93957 * (A - Z) - 0.001 * bindingEnergy(A, Z);
  };

  void setExitonConfiguration(const G4ExitonConfiguration& config) { 

    theExitonConfiguration = config;
  };

  G4double getA() const { 

    return A; 
  };

  G4double getZ() const { 

    return Z; 
  };

  G4double getExitationEnergy() const { 

    return exitationEnergy; 
  };

  G4double getExitationEnergyInGeV() const { 

    return 0.001 * exitationEnergy; 
  };

  G4ExitonConfiguration getExitonConfiguration() const {

    return theExitonConfiguration;
  };

  G4double getMass() const { 

    return nucleiMass; 
  };

  G4double getKineticEnergy() const { 

    return momentum[0] - nucleiMass; 
  };

  G4double getNucleiMass(G4double a, 
			 G4double z) const {

    return 0.93827 * z + 0.93957 * (a - z) - 0.001 * bindingEnergy(a, z);
  };

  virtual void printParticle() const {

    G4cout << " A " << A << " Z " << Z << " mass " << nucleiMass << 
      " Eex (MeV) " << exitationEnergy << G4endl;

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






