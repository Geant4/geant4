// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
#include "G4NeutronHPSCFissionFS.hh"

  void G4NeutronHPSCFissionFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
  {
    G4String aString = "/SC/";
    G4NeutronHPFissionBaseFS::Init(A, Z, dirName, aString);
  }
  
  G4DynamicParticleVector * G4NeutronHPSCFissionFS::ApplyYourself(G4int NNeutrons)
  {  
    G4DynamicParticleVector * aResult;
//    G4cout <<"G4NeutronHPSCFissionFS::ApplyYourself +"<<endl;
    aResult = G4NeutronHPFissionBaseFS::ApplyYourself(NNeutrons);    
    return aResult;
  }
