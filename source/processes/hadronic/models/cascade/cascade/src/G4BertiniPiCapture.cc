#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4Nucleus.hh"

#include "G4BertiniPiCapture.hh"

G4BertiniPiCapture::G4BertiniPiCapture() {
  ;
}

G4BertiniPiCapture::~G4BertiniPiCapture(){
  ;
}

/*
// pi- capture
    if (ifirst == 1) {
      G4int m;
      for(m = 0; m < mxmat; ++m) {
	lmx = nel[m];
	for(l = 0; l < lmx; ++l) pcap[l][m] = den[l][m] * zz[l][m];
	pcaps = 0.0;
	for(l = 0; l < lmx; ++l) pcaps += pcap[l][m];
	for(l = 0; l < lmx; ++l) pcap[l][m] /= pcaps;
      }
      ifirst = 0;
    }
    lmx = nel[mat];
    G4bool flag = true;
    for(l = 0; l < lmx; ++l) {
      r0 = r;
      r -= pcap[l][mat];
      if (r <= 0.0) {
	lelem = l;
	flag = false;
	break;
      }
    }
    if (flag) G4cerr << "cascad1" << G4endl;
    dkwt = 1.0;
    if (lelem == 0 || nhist == nf) {
      if (verboseLevel > 1) {
	G4cout << " cascad "          << G4endl;
	G4cout << " lelem  " << lelem << G4endl;
	G4cout << " mat    " << mat   << G4endl;
	G4cout << " nhist  " << nhist << G4endl;
      }
    }

*/





