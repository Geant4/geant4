// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHP4NInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Alpha.hh"

  G4NeutronHP4NInelasticFS::G4NeutronHP4NInelasticFS(){}
  G4NeutronHP4NInelasticFS::~G4NeutronHP4NInelasticFS(){}
G4ParticleChange * G4NeutronHP4NInelasticFS::ApplyYourself(const G4Track & theTrack)
{
// these are the particle types in the final state

  G4ParticleDefinition * theDefs[4];
  theDefs[0] = G4Neutron::Neutron();
  theDefs[1] = G4Neutron::Neutron();
  theDefs[2] = G4Neutron::Neutron();
  theDefs[3] = G4Neutron::Neutron();
  
// fill the final state  
  G4NeutronHPInelasticBaseFS::BaseApply(theTrack, theDefs, 4);
  
// return the result
   return &theResult;
}

void G4NeutronHP4NInelasticFS::
Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticBaseFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-3;
   G4double ResidualZ = Z;
   G4NeutronHPInelasticBaseFS::InitGammas(ResidualA, ResidualZ);
}
