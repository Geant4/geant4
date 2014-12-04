// -------------------------------------------------------------
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFPrimaryGeneratorAction.cc
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Implementation of the FFPrimaryGeneratorAction class
//!
//! \details    Generates 4.5 MeV neutrons from the solid volume defined in
//!             FFDetectorConstruction as "NeturonSource" and shoots them off
//!             in a isotropically sampled direction
//!
//  ================ End Documentation Comments ================
//
//  Modified:
//
//  23-06-14                                              BWendt
//  Implemented "GetNeutronSourceCenter()" and added code to correctly sample
//  the neutron starting location
//
// -------------------------------------------------------------

#include "globals.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Neutron.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"

#include "FFPrimaryGeneratorAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFPrimaryGeneratorAction::
FFPrimaryGeneratorAction()
:   G4VUserPrimaryGeneratorAction(),
#ifndef NDEBUG
    fEventNumber(0),
#endif // NDEBUG
    fH2OPhysical(NULL),
    fNeutronPhysical(NULL),
    fNeutronSolid(NULL),
    fParticleGun(new G4ParticleGun(1)),
    fTankPhysical(NULL)
{
    fParticleGun->SetParticleDefinition(G4Neutron::Definition());
    fParticleGun->SetParticleEnergy(4.5 * MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FFPrimaryGeneratorAction::
GeneratePrimaries(G4Event* event)
{
#ifndef NDEBUG
    G4cout << "Shooting event " << ++fEventNumber << G4endl;
#endif // NDEBUG
    const G4ThreeVector sourceCenter = GetNeutronSourceCenter();
    
    // Sample the neutron source location
    const G4double radius = fNeutronSolid->GetOuterRadius();
    const G4double z = fNeutronSolid->GetZHalfLength() * 2;
    G4ThreeVector randomLocation;
    randomLocation.setRThetaPhi(radius * std::sqrt(G4UniformRand()),
                                G4UniformRand() * 180 * deg,
                                0);
    randomLocation.setZ(z * (G4UniformRand() - 0.5));
    G4ThreeVector location(randomLocation.x() + sourceCenter.x(),
                           randomLocation.y() + sourceCenter.y(),
                           randomLocation.z() + sourceCenter.z());
#ifndef NDEBUG
    G4cout << "Emission Location: r: " << location << G4endl;
#endif // NDEBUG
    
    // Sample the neutron emission direction
    G4ThreeVector direction;
    direction.setRThetaPhi(1.0,
                           std::acos(G4UniformRand() * 2 - 1),
                           (G4UniformRand() * 2 - 1) * 180 * deg);
#ifndef NDEBUG
    G4cout << "Emission Direction: r: " << direction << G4endl;
#endif // NDEBUG
    
    // Load the event
    fParticleGun->SetParticlePosition(location);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector FFPrimaryGeneratorAction::
GetNeutronSourceCenter(void)
{
    // Get the dimensions of the neutron source
    if(fNeutronSolid == NULL)
    {
        G4LogicalVolume* temp
            = G4LogicalVolumeStore::GetInstance()->GetVolume("NeutronSource");
        if(temp != NULL)
        {
            fNeutronSolid = dynamic_cast< G4Tubs* >(temp->GetSolid());
        }
        
        if(fNeutronSolid == NULL)
        {
            G4Exception("FFPrimaryGeneratorAction::"
                            "GeneratePrimaries(G4Event*)",
                        "Neutron source solid volume not found",
                        EventMustBeAborted,
                        "This run will be aborted");
        }
    }
    
    // Get the position of the neutron source within the water
    if(fNeutronPhysical == NULL)
    {
        fNeutronPhysical = G4PhysicalVolumeStore::GetInstance()
            ->GetVolume("NeutronSource");
    }
    
    if(fNeutronPhysical == NULL)
    {
        G4Exception("FFPrimaryGeneratorAction::GetNeutronSourceCenter(void)",
                    "Neutron source physical volume not found",
                    EventMustBeAborted,
                    "This run will be aborted");
    }
    
    // Get the position of the water within the tank
    if(fH2OPhysical == NULL)
    {
        fH2OPhysical = G4PhysicalVolumeStore::GetInstance()
            ->GetVolume("Tank_H2O");
    
        if(fH2OPhysical == NULL)
        {
            G4Exception("FFPrimaryGeneratorAction::"
                            "GetNeutronSourceCenter(void)",
                        "Tank H2O physical volume not found",
                        EventMustBeAborted,
                        "This run will be aborted");
        }
    }
    
    // Get the position of the tank within the world
    if(fTankPhysical == NULL)
    {
        fTankPhysical = G4PhysicalVolumeStore::GetInstance()
            ->GetVolume("Tank_Wall");
    
        if(fTankPhysical == NULL)
        {
            G4Exception("FFPrimaryGeneratorAction::"
                            "GetNeutronSourceCenter(void)",
                        "Tank physical volume not found",
                        EventMustBeAborted,
                        "This run will be aborted");
        }
    }
    
    return fNeutronPhysical->GetTranslation()
               + fH2OPhysical->GetTranslation()
               + fTankPhysical->GetTranslation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFPrimaryGeneratorAction::
~FFPrimaryGeneratorAction()
{
    delete fParticleGun;
}


