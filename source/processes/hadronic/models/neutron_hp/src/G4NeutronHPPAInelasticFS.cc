// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPPAInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Alpha.hh"
#include "G4Proton.hh"

G4ParticleChange * G4NeutronHPPAInelasticFS::ApplyYourself(const G4Track & theTrack)
{
// these are the particle types in the final state

  G4ParticleDefinition * theDefs[2];
  theDefs[0] = G4Proton::Proton();
  theDefs[1] = G4Alpha::Alpha();
  
// fill the final state  
  G4NeutronHPInelasticBaseFS::BaseApply(theTrack, theDefs, 2);
  
// return the result
   return &theResult;
}

void G4NeutronHPPAInelasticFS::
Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticBaseFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-4;
   G4double ResidualZ = Z-3;
   G4NeutronHPInelasticBaseFS::InitGammas(ResidualA, ResidualZ);
}
