// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPDInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Deuteron.hh"

void G4NeutronHPDInelasticFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticCompFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-1;
   G4double ResidualZ = Z-1;
   G4NeutronHPInelasticCompFS::InitGammas(ResidualA, ResidualZ);
}

G4ParticleChange * G4NeutronHPDInelasticFS::ApplyYourself(const G4Track & theTrack)
{

// do the final state
    G4NeutronHPInelasticCompFS::CompositeApply(theTrack, G4Deuteron::Deuteron());
             
// return the result
    return &theResult;
}
