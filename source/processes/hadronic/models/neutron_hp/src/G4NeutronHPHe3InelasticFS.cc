// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPHe3InelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4He3.hh"

void G4NeutronHPHe3InelasticFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticCompFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A-2;
   G4double ResidualZ = Z-2;
   G4NeutronHPInelasticCompFS::InitGammas(ResidualA, ResidualZ);
}

G4ParticleChange * G4NeutronHPHe3InelasticFS::ApplyYourself(const G4Track & theTrack)
{

// do the final state
    G4NeutronHPInelasticCompFS::CompositeApply(theTrack, G4He3::He3());
             
// return the result
    return &theResult;
}
