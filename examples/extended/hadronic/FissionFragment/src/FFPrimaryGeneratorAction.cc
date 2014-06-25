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
    H2OPhysical(NULL),
    neutronPhysical(NULL),
    neutronSolid(NULL),
    particleGun(new G4ParticleGun(1)),
    tankPhysical(NULL)
{
    particleGun->SetParticleDefinition(G4Neutron::Definition());
    particleGun->SetParticleEnergy(4.5 * MeV);
    //particleGun->SetParticleEnergy(10.5 * MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FFPrimaryGeneratorAction::
GeneratePrimaries(G4Event* event)
{
    const G4ThreeVector sourceCenter = GetNeutronSourceCenter();
    
    // Sample the neutron source location
    const G4double radius = neutronSolid->GetOuterRadius();
    const G4double z = neutronSolid->GetZHalfLength() * 2;
    G4ThreeVector randomLocation;
    randomLocation.setRThetaPhi(radius * sqrt(G4UniformRand()),
                                G4UniformRand() * 180 * deg,
                                0);
    randomLocation.setZ(z * (G4UniformRand() - 0.5));
    G4ThreeVector location(randomLocation.x() + sourceCenter.x(),
                           randomLocation.y() + sourceCenter.y(),
                           randomLocation.z() + sourceCenter.z());
    //G4cout << "Emission Location: r: " << location << G4endl;
    
    // Sample the neutron emission direction
    G4ThreeVector direction;
    direction.setRThetaPhi(1.0,
                           std::acos(G4UniformRand() * 2 - 1),
                           (G4UniformRand() * 2 - 1) * 180 * deg);
    //G4cout << "Emission Direction: r: " << direction << G4endl;
    
    // Load the event
    particleGun->SetParticlePosition(location);
    particleGun->SetParticleMomentumDirection(direction);
    particleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector FFPrimaryGeneratorAction::
GetNeutronSourceCenter(void)
{
    // Get the dimensions of the neutron source
    if(neutronSolid == NULL)
    {
        G4LogicalVolume* temp
            = G4LogicalVolumeStore::GetInstance()->GetVolume("NeutronSource");
        if(temp != NULL)
        {
            neutronSolid = dynamic_cast< G4Tubs* >(temp->GetSolid());
        }
        
        if(neutronSolid == NULL)
        {
            G4Exception("FFPrimaryGeneratorAction::"
                            "GeneratePrimaries(G4Event*)",
                        "Neutron source solid volume not found",
                        EventMustBeAborted,
                        "This run will be aborted");
        }
    }
    
    // Get the position of the neutron source within the water
    if(neutronPhysical == NULL)
    {
        neutronPhysical = G4PhysicalVolumeStore::GetInstance()
            ->GetVolume("NeutronSource");
    }
    
    if(neutronPhysical == NULL)
    {
        G4Exception("FFPrimaryGeneratorAction::GetNeutronSourceCenter(void)",
                    "Neutron source physical volume not found",
                    EventMustBeAborted,
                    "This run will be aborted");
    }
    
    // Get the position of the water within the tank
    if(H2OPhysical == NULL)
    {
        H2OPhysical = G4PhysicalVolumeStore::GetInstance()
            ->GetVolume("Tank_H2O");
    
        if(H2OPhysical == NULL)
        {
            G4Exception("FFPrimaryGeneratorAction::"
                            "GetNeutronSourceCenter(void)",
                        "Tank H2O physical volume not found",
                        EventMustBeAborted,
                        "This run will be aborted");
        }
    }
    
    // Get the position of the tank within the world
    if(tankPhysical == NULL)
    {
        tankPhysical = G4PhysicalVolumeStore::GetInstance()
            ->GetVolume("Tank_Wall");
    
        if(tankPhysical == NULL)
        {
            G4Exception("FFPrimaryGeneratorAction::"
                            "GetNeutronSourceCenter(void)",
                        "Tank physical volume not found",
                        EventMustBeAborted,
                        "This run will be aborted");
        }
    }
    
    return neutronPhysical->GetTranslation()
               + H2OPhysical->GetTranslation()
               + tankPhysical->GetTranslation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFPrimaryGeneratorAction::
~FFPrimaryGeneratorAction()
{
    delete particleGun;
}


