// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPDAInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"

  G4NeutronHPDAInelasticFS::G4NeutronHPDAInelasticFS(){}
  G4NeutronHPDAInelasticFS::~G4NeutronHPDAInelasticFS(){}
G4ParticleChange * G4NeutronHPDAInelasticFS::ApplyYourself(const G4Track & theTrack)
{
// these are the particle types in the final state

  G4ParticleDefinition * theDefs[2];
  theDefs[0] = G4Deuteron::Deuteron();
  theDefs[1] = G4Alpha::Alpha();
  
// fill the final state  
  G4NeutronHPInelasticBaseFS::BaseApply(theTrack, theDefs, 2);
  
// return the result
   return &theResult;
}

void G4NeutronHPDAInelasticFS::
Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticBaseFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-5;
   G4double ResidualZ = Z-3;
   G4NeutronHPInelasticBaseFS::InitGammas(ResidualA, ResidualZ);
}
