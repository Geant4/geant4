#include "G4ImportanceAlgorithm.hh"
#include "Randomize.hh"
#include "g4std/strstream"

G4ImportanceAlgorithm::G4ImportanceAlgorithm(): fWorned(false) {}
G4ImportanceAlgorithm::~G4ImportanceAlgorithm(){
  if(fWorned) {
    G4cout << G4endl;
    Warning("~G4ImportanceAlgorithm: ipre_over_ipost ! in [0.25, 4] seen");
    G4cout << G4endl;
  }
}

G4Nsplit_Weight
G4ImportanceAlgorithm::
Calculate(G4double ipre_over_ipost, G4double init_w) const{

  if ((ipre_over_ipost<0.25 || ipre_over_ipost> 4) && !fWorned) {
    ostrstream os;
    os << "Calculate: ipre_over_ipost ! in [0.25, 4]: ipre_over_ipost = "
       << ipre_over_ipost << '\0' << G4endl;
    Warning(os.str());
    fWorned = true;
    if (ipre_over_ipost<=0) {
      Error("Calculate: ipre_over_ipost<=0");
    }
  }
  if (init_w<=0.) {
    Error("Calculate:  iniitweight<= 0. found");
  }

  // default geometrical splitting 
  // in integer mode 
  // for ipre_over_ipost <= 1
  G4double inv = 1./ipre_over_ipost;
  G4Nsplit_Weight nw( (int)inv, init_w * ipre_over_ipost);
  
  // geometrical splitting for double mode
  if (ipre_over_ipost<1) {
    if ( G4double(nw.fN) != inv) {
      // double mode
      // probability p for splitting into n+1 tracks
      G4double p = inv - nw.fN;
      // get a random number out of [0,1)
      G4double r = G4UniformRand();
      if (r<p) {
	nw.fN++;
      } 
    }  
  }
  // ipre_over_ipost > 1
  //  russian roulett
  else {
    // probabiity for killing track
    G4double p = 1-inv;
    // get a random number out of [0,1)
    G4double r = G4UniformRand();
    if (r<p) {
      // kill track
      nw.fN = 0;
      nw.fW = 0;
    }
    else {
      nw.fN = 1;     
    }
  }
  

  return nw;
}
