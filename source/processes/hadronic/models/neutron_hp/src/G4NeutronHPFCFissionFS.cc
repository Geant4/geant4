// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPFCFissionFS.hh"

  void G4NeutronHPFCFissionFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
  {
    G4String aString = "/FC/";
    G4NeutronHPFissionBaseFS::Init(A, Z, dirName, aString);
  }
  
  G4DynamicParticleVector * G4NeutronHPFCFissionFS::ApplyYourself(G4int nNeutrons)
  {  
    G4DynamicParticleVector * aResult;
//    G4cout <<"G4NeutronHPFCFissionFS::ApplyYourself +"<<G4endl;
    aResult = G4NeutronHPFissionBaseFS::ApplyYourself(nNeutrons);    
    return aResult;
  }
