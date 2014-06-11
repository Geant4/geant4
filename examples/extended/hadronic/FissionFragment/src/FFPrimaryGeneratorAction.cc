// -------------------------------------------------------------
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFPrimaryGeneratorAction.hh
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Implementation of the FFPrimaryGeneratorAction class
//!
//! \details    Generates 4.5 MeV neutrons from the solid volume defined in
//!             FFDetectorConstruction as "NeturonSource" and shoots them off
//!             in a uniformly sampled direction over the surface of a sphere
//!
//  ================ End Documentation Comments ================
//
//  Modified: 
//
// -------------------------------------------------------------

#include "globals.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Neutron.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"

#include "FFPrimaryGeneratorAction.hh"


FFPrimaryGeneratorAction::
FFPrimaryGeneratorAction()
:   G4VUserPrimaryGeneratorAction(),
    neutronSource(NULL),
    particleGun(new G4ParticleGun(1))
{
    particleGun->SetParticleDefinition(G4Neutron::Definition());
    particleGun->SetParticleEnergy(4.5 * MeV);
}

void FFPrimaryGeneratorAction::
GeneratePrimaries(G4Event* event)
{
    if(neutronSource == NULL)
    {
        G4LogicalVolume* temp
            = G4LogicalVolumeStore::GetInstance()->GetVolume("NeutronSource");
        if(temp != NULL)
        {
            neutronSource = dynamic_cast< G4Tubs* >(temp->GetSolid());
        }
    }
    
    if(neutronSource == NULL)
    {
        G4Exception("FFPrimaryGeneratorAction::GeneratePrimaries(G4Event*)",
                    "Neutron source geometry not found",
                    EventMustBeAborted,
                    "This run will be aborted");
    }
    
    // Sample the neutron source location
    const G4double radius = neutronSource->GetOuterRadius();
    const G4double z = neutronSource->GetZHalfLength() * 2;
    G4ThreeVector location;
    location.setRThetaPhi(radius * sqrt(G4UniformRand()),
                          G4UniformRand() * 180 * deg,
                          0);
    location.setZ(z * (G4UniformRand() - 0.5));
    
    // Sample the neutron emission direction
    G4ThreeVector direction;
    direction.setRThetaPhi(1.0,
                           std::acos(G4UniformRand() * 2 - 1),
                           (G4UniformRand() * 2 - 1) * 180 * deg);
    
    // Load the event
    particleGun->SetParticlePosition(location);
    particleGun->SetParticleMomentumDirection(direction);
    particleGun->GeneratePrimaryVertex(event);
}

FFPrimaryGeneratorAction::
~FFPrimaryGeneratorAction()
{
    delete particleGun;
}


