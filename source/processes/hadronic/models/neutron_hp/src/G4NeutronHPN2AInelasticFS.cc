// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPN2AInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Alpha.hh"

G4ParticleChange * G4NeutronHPN2AInelasticFS::ApplyYourself(const G4Track & theTrack)
{
// these are the particle types in the final state

  G4ParticleDefinition * theDefs[3];
  theDefs[0] = G4Neutron::Neutron();
  theDefs[1] = G4Alpha::Alpha();
  theDefs[2] = G4Alpha::Alpha();
  
// fill the final state  
  G4NeutronHPInelasticBaseFS::BaseApply(theTrack, theDefs, 3);
  
// return the result
   return &theResult;
}

void G4NeutronHPN2AInelasticFS::
Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticBaseFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-8;
   G4double ResidualZ = Z-4;
   G4NeutronHPInelasticBaseFS::InitGammas(ResidualA, ResidualZ);
}

