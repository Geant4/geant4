// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPNInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/strstream"

void G4NeutronHPNInelasticFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
   G4NeutronHPInelasticCompFS::Init(A, Z, dirName, aFSType);
   G4double ResidualA = A;
   G4double ResidualZ = Z;
   G4NeutronHPInelasticCompFS::InitGammas(ResidualA, ResidualZ);
}

G4ParticleChange * G4NeutronHPNInelasticFS::ApplyYourself(const G4Track & theTrack)
{

// do the final state
    G4NeutronHPInelasticCompFS::CompositeApply(theTrack, G4Neutron::Neutron());
             
// return the result
    return &theResult;
}
