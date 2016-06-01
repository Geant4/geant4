// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPD2AInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"

G4ParticleChange * G4NeutronHPD2AInelasticFS::ApplyYourself(const G4Track & theTrack)
{
// these are the particle types in the final state

  G4ParticleDefinition * theDefs[3];
  theDefs[0] = G4Deuteron::Deuteron();
  theDefs[1] = G4Alpha::Alpha();
  theDefs[2] = G4Alpha::Alpha();
  
// fill the final state  
  G4NeutronHPInelasticBaseFS::BaseApply(theTrack, theDefs, 3);
  
// return the result
   return &theResult;
}

void G4NeutronHPD2AInelasticFS::
Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticBaseFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-9;
   G4double ResidualZ = Z-5;
   G4NeutronHPInelasticBaseFS::InitGammas(ResidualA, ResidualZ);
}
