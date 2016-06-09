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
#ifndef G4EXITON_CONFIGURATION_HH
#define G4EXITON_CONFIGURATION_HH

class G4ExitonConfiguration {

public:

  G4ExitonConfiguration() {
    protonQuasiParticles = 0.0;
    neutronQuasiParticles = 0.0;
    protonHoles = 0.0;
    neutronHoles = 0.0;
  };

  G4ExitonConfiguration(G4double qpp, 
			G4double qnp, 
			G4double qph, 
			G4double qnh) 
    : protonQuasiParticles(qpp), 
    neutronQuasiParticles(qnp), 
    protonHoles(qph), 
    neutronHoles(qnh) {
  };
 
  void incrementQP(G4int ip) {
    if(ip < 3) {
      if(ip == 1) {
	protonQuasiParticles += 1.0;
      }
      else if(ip == 2) {
	neutronQuasiParticles += 1.0;
      };
    };
  };

  void incrementHoles(G4int ip) {
    if(ip < 3) {
      if(ip == 1) {
	protonHoles += 1.0;
      }
      else if(ip == 2) {
	neutronHoles += 1.0;
      };
    };
  };

  void print() const {
    G4cout << " Exiton configuration " << G4endl
	   << " proton particles " << protonQuasiParticles << " holes " 
	   << protonHoles << G4endl
	   << " neutron particles " << neutronQuasiParticles << " holes " 
	   << neutronHoles << G4endl;
  };
     
  G4double protonQuasiParticles;
  G4double neutronQuasiParticles;
  G4double protonHoles;
  G4double neutronHoles;

};        

#endif // G4EXITON_CONFIGURATION_HH 
