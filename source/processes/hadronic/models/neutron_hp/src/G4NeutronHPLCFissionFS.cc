// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPLCFissionFS.hh"

  G4NeutronHPLCFissionFS::G4NeutronHPLCFissionFS(){ hasXsec = false; }
  G4NeutronHPLCFissionFS::~G4NeutronHPLCFissionFS(){}
  void G4NeutronHPLCFissionFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
  {
    G4String aString = "/LC/";
    G4NeutronHPFissionBaseFS::Init(A, Z, dirName, aString);
  }
  
  G4DynamicParticleVector * G4NeutronHPLCFissionFS::ApplyYourself(G4int NNeutrons)
  {  
    G4DynamicParticleVector * aResult;
//    G4cout <<"G4NeutronHPLCFissionFS::ApplyYourself +"<<endl;
    aResult = G4NeutronHPFissionBaseFS::ApplyYourself(NNeutrons);    
    return aResult;
  }
