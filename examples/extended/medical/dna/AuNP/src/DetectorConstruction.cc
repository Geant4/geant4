//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file medical/dna/range/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 78723 2014-01-20 10:32:17Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4ProductionCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fNPRadius(0),
   fAbsRadius(0),
   fNPMaterial(0),
   fAbsMaterial(0),
   pWorld(0),
   fNP(0),
   fAbs(0),
   fDetectorMessenger(0),
   fRegion(0)
{
  //default tracking cut  
  fTrackingCut = 10.0*eV;
  
  // default parameter values
  fNPRadius    = 50  *nm;
  fAbsRadius   = 100000 *nm+fNPRadius;

  fNreplicaR   = 1000;
  fNreplicaAzm = 360;
  
  DefineMaterials();
  SetNPMaterial ("G4_Au");
  SetAbsMaterial("G4_WATER");

  // create commands for interactive definition of the detector  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* man = G4NistManager::Instance();
  
  man->FindOrBuildMaterial("G4_Au");
  man->FindOrBuildMaterial("G4_WATER");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  const G4bool check_overlap = false;
  G4GeometryManager::GetInstance()    ->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance() ->Clean();
  G4SolidStore::GetInstance()         ->Clean();

  //World  ====================================================================
  G4Material* Vacuum =
          new G4Material("Galactic", 1., 1.01*g/mole,universe_mean_density,
                         kStateGas, 2.73*kelvin, 3.e-18*pascal);

  G4Sphere         *sWorld = new G4Sphere("World",0., fAbsRadius*1.5, 0.,
                                          2*pi*rad, 0., pi);
  G4LogicalVolume  *lWorld = new G4LogicalVolume(sWorld,Vacuum,"World");
                    pWorld = new G4PVPlacement(0,G4ThreeVector(),lWorld,
                                               "World",0,false,0);

                   
  // Spherical NanoParticle ===================================================
  G4Sphere* 
  sNP        = new G4Sphere("NanoParticle",                      
                            0., fNPRadius, 0., 2*pi*rad, 0., pi);

  fLogicalNP = new G4LogicalVolume(sNP,               
                                   fNPMaterial,       
                                   "NanoParticle");   
                                   
  fNP        = new G4PVPlacement(0,                   
                                 G4ThreeVector(),     
                                 fLogicalNP,          
                                 "NanoParticle",      
                                 lWorld,              
                                 false,               
                                 0,                   
                                 check_overlap);      

  // Sampling Plane ===========================================================
  
  // Spherical water absorber -------------------------------------------------
  G4Sphere* 
  sAbsorber0 = new G4Sphere("Absorber0",                          
                            0., fAbsRadius*1.2, 0., 2*pi*rad, 0., pi);
  G4VSolid*
  sAbsorber = new G4SubtractionSolid("Absorber",
                                      sAbsorber0,sNP,0,G4ThreeVector());

  fLogicalAbs= new G4LogicalVolume(sAbsorber,     
                                   fAbsMaterial,  
                                   "Absorber");   
  fAbs        = new G4PVPlacement(0,               
                                  G4ThreeVector(), 
                                  fLogicalAbs,     
                                  "Absorber",      
                                  lWorld,          
                                  false,           
                                  0,               
                                  check_overlap);  
  

  //Set Visualization Setting#################################################
  G4VisAttributes *visWorld = new G4VisAttributes(true,G4Colour::Gray());
  G4VisAttributes *visNP    = new G4VisAttributes(true,G4Colour::Yellow());
  G4VisAttributes *visAbs   = new G4VisAttributes(true,G4Colour::Cyan());

  lWorld     ->SetVisAttributes(visWorld);
  fLogicalNP ->SetVisAttributes(visNP);
  fLogicalAbs->SetVisAttributes(visAbs);

  PrintParameters();
    
  G4double maxStep = fNPRadius/2.;
  fLogicalNP->SetUserLimits(new G4UserLimits(maxStep,DBL_MAX,DBL_MAX,
                            fTrackingCut));
                            

  fRegion = new G4Region("NP");

  G4ProductionCuts *cuts = new G4ProductionCuts();
  G4double defCut = 0.1*nm;
  cuts->SetProductionCut(defCut,"gamma");
  cuts->SetProductionCut(defCut,"e-");
  cuts->SetProductionCut(defCut,"e+");
  cuts->SetProductionCut(defCut,"proton");
  fRegion->SetProductionCuts(cuts);
  fRegion->AddRootLogicalVolume(fLogicalNP);
  
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() const
{
  G4cout << "\n=========================================================\n";
  G4cout << "---> The tracking cut to all particles is set to " 
         << G4BestUnit(fTrackingCut,"Energy") << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Nano-Particle is a sphere of " 
         << G4BestUnit(fNPRadius,"Length") << " radius of "
         << fNPMaterial->GetName() << " made of"
         << "\n \n" << fNPMaterial << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is a sphere of " 
         << G4BestUnit(fAbsRadius,"Length") << " radius of "
         << fAbsMaterial->GetName() << " made of"
         << "\n \n" << fAbsMaterial << G4endl;
  G4cout << "\n=========================================================\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTrackingCut(G4double value)
{
  fTrackingCut = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsRadius(G4double value)
{
  fAbsRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}

void DetectorConstruction::SetNPRadius(G4double value)
{
  fNPRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
void DetectorConstruction::SetAbsMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) fAbsMaterial = pttoMaterial;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();  
}

void DetectorConstruction::SetNPMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) fNPMaterial = pttoMaterial;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();  
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
void DetectorConstruction::SetNReplicaR(G4int value)
{
  fNreplicaR = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}
void DetectorConstruction::SetNReplicaAzm(G4int value)
{
  fNreplicaAzm = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}
