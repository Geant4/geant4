// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
#include "G4NeutronHPTCFissionFS.hh"

  G4NeutronHPTCFissionFS::G4NeutronHPTCFissionFS(){ hasXsec = false; }
  G4NeutronHPTCFissionFS::~G4NeutronHPTCFissionFS(){}

  void G4NeutronHPTCFissionFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType0)
  {
    G4String aString = "/TC/";
    G4NeutronHPFissionBaseFS::Init(A, Z, dirName, aString);
  }
  
  G4DynamicParticleVector * G4NeutronHPTCFissionFS::ApplyYourself(G4int NNeutrons)
  {  
    G4DynamicParticleVector * aResult;
//    G4cout <<"G4NeutronHPTCFissionFS::ApplyYourself +"<<endl;
    aResult = G4NeutronHPFissionBaseFS::ApplyYourself(NNeutrons);    
    return aResult;
  }
