// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPNHe3InelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4He3.hh"

G4ParticleChange * G4NeutronHPNHe3InelasticFS::ApplyYourself(const G4Track & theTrack)
{
// these are the particle types in the final state

  G4ParticleDefinition * theDefs[2];
  theDefs[0] = G4Neutron::Neutron();
  theDefs[1] = G4He3::He3();
  
// fill the final state  
  G4NeutronHPInelasticBaseFS::BaseApply(theTrack, theDefs, 2);
  
// return the result
   return &theResult;
}

void G4NeutronHPNHe3InelasticFS::
Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticBaseFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-3;
   G4double ResidualZ = Z-2;
   G4NeutronHPInelasticBaseFS::InitGammas(ResidualA, ResidualZ);
}
