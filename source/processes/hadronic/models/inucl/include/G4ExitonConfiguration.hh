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
