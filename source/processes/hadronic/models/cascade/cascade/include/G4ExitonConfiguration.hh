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
