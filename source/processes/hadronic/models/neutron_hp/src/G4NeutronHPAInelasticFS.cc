// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPAInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Alpha.hh"

void G4NeutronHPAInelasticFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticCompFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-3;
   G4double ResidualZ = Z-2;
   G4NeutronHPInelasticCompFS::InitGammas(ResidualA, ResidualZ);
}

G4ParticleChange * G4NeutronHPAInelasticFS::ApplyYourself(const G4Track & theTrack)
{

// do the final state
   G4NeutronHPInelasticCompFS::CompositeApply(theTrack, G4Alpha::Alpha());
             
// return the result
   return &theResult;
}
