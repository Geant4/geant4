// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPPInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"

void G4NeutronHPPInelasticFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticCompFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A;
   G4double ResidualZ = Z-1;
   G4NeutronHPInelasticCompFS::InitGammas(ResidualA, ResidualZ);
}

G4ParticleChange * G4NeutronHPPInelasticFS::ApplyYourself(const G4Track & theTrack)
{

// do the final state
    G4NeutronHPInelasticCompFS::CompositeApply(theTrack, G4Proton::Proton());
             
// return the result
    return &theResult;
}

