#ifndef EXITON_CONFIGURATION_H
#define EXITON_CONFIGURATION_H

class ExitonConfiguration {

public:

ExitonConfiguration() {
  protonQuasiParticles = neutronQuasiParticles = protonHoles = neutronHoles 
   = 0.;
};

ExitonConfiguration(double qpp, double qnp, double qph, double qnh) :
  protonQuasiParticles(qpp), neutronQuasiParticles(qnp), 
  protonHoles(qph), neutronHoles(qnh) {};
 
void incrementQP(int ip) {
  if(ip < 3) {
    if(ip == 1) {
      protonQuasiParticles += 1.;
    }
     else if(ip == 2) {
      neutronQuasiParticles += 1.;
    };
  };
};

void incrementHoles(int ip) {
  if(ip < 3) {
    if(ip == 1) {
      protonHoles += 1.;
    }
     else if(ip == 2) {
      neutronHoles += 1.;
    };
  };
};

void print() const {
  cout << " Exiton configuration " << endl
       << " proton particles " << protonQuasiParticles << " holes " 
       << protonHoles << endl
       << " neutron particles " << neutronQuasiParticles << " holes " 
       << neutronHoles << endl;
};
     
double protonQuasiParticles;
double neutronQuasiParticles;
double protonHoles;
double neutronHoles;

};        

#endif // EXITON_CONFIGURATION_H 
