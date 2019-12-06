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
/// \file electromagnetic/TestEm2/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4GlobalMagFieldMessenger.hh"

#include "G4NistManager.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4AutoDelete.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 :G4VUserDetectorConstruction(),
  fNLtot(40),fNRtot(50),fDLradl(0.5),fDRradl(0.1),
  fDLlength(0.),fDRlength(0.),fMaterial(nullptr),
  fEcalLength(0.),fEcalRadius(0.),
  fSolidEcal(nullptr),fLogicEcal(nullptr),fPhysiEcal(nullptr)
{
  DefineMaterials();
  SetMaterial("G4_PbWO4");
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //
  // define few Elements by hand
  //
  G4double a, z;
    
  G4Element* H  = new G4Element("Hydrogen",  "H", z= 1., a=   1.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,  "O", z= 8., a=  16.00*g/mole);
  G4Element* Ge = new G4Element("Germanium", "Ge",z=32., a=  72.59*g/mole);
  G4Element* Bi = new G4Element("Bismuth",   "Bi",z=83., a= 208.98*g/mole);

  //
  // define materials
  //
  G4double density;
  G4int ncomponents, natoms;

  // water with ionisation potential 78 eV
  G4Material* H2O = 
  new G4Material("Water", density= 1.00*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  // pure materails
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  new G4Material("Aluminium",   z=13., a= 26.98*g/mole, density= 2.7*g/cm3);
  new G4Material("Iron",        z=26., a= 55.85*g/mole, density= 7.87*g/cm3);  
  new G4Material("Copper",      z=29., a= 63.55*g/mole, density= 8.960*g/cm3); 
  new G4Material("Tungsten",    z=74., a=183.84*g/mole, density=19.35*g/cm3); 
  new G4Material("Lead",        z=82., a=207.19*g/mole, density=11.35*g/cm3);  
  new G4Material("Uranium"    , z=92., a=238.03*g/mole, density= 18.95*g/cm3);

  // compound material
  G4Material* BGO = 
  new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
  BGO->AddElement(O , natoms=12);
  BGO->AddElement(Ge, natoms= 3);
  BGO->AddElement(Bi, natoms= 4);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateParameters()
{
  G4double Radl = fMaterial->GetRadlen();
  fDLlength = fDLradl*Radl; fDRlength = fDRradl*Radl;
  fEcalLength = fNLtot*fDLlength;  fEcalRadius = fNRtot*fDRlength;
  if(fSolidEcal) {
    fSolidEcal->SetOuterRadius(fEcalRadius);
    fSolidEcal->SetZHalfLength(0.5*fEcalLength);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  UpdateParameters();
  //
  // Ecal
  //
  if(!fPhysiEcal) { 
    fSolidEcal = new G4Tubs("Ecal",0.,fEcalRadius,0.5*fEcalLength,0.,360*deg);
    fLogicEcal = new G4LogicalVolume( fSolidEcal,fMaterial,"Ecal",0,0,0);
    fPhysiEcal = new G4PVPlacement(0,G4ThreeVector(),
                                   fLogicEcal,"Ecal",0,false,0);
  }
  G4cout << "Absorber is " << G4BestUnit(fEcalLength,"Length")
         << " of " << fMaterial->GetName()
         << "  R= " << fEcalRadius/cm << " cm" << G4endl;
  //
  //always return the physical World
  //
  return fPhysiEcal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if(pttoMaterial &&  fMaterial != pttoMaterial) {
    fMaterial = pttoMaterial;
    if(fLogicEcal) { fLogicEcal->SetMaterial(fMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetLBining(G4ThreeVector Value)
{
  fNLtot = (G4int)Value(0);
  if (fNLtot > kMaxBin) {
    G4cout << "\n ---> warning from SetLBining: "
           << fNLtot << " truncated to " << kMaxBin << G4endl;
    fNLtot = kMaxBin;
  }  
  fDLradl = Value(1);
  UpdateParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRBining(G4ThreeVector Value)
{
  fNRtot = (G4int)Value(0);
  if (fNRtot > kMaxBin) {
    G4cout << "\n ---> warning from SetRBining: "
           << fNRtot << " truncated to " << kMaxBin << G4endl;
    fNRtot = kMaxBin;
  }    
  fDRradl = Value(1);
  UpdateParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  if ( fFieldMessenger.Get() == nullptr ) {
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    G4GlobalMagFieldMessenger* msg =
      new G4GlobalMagFieldMessenger(fieldValue);
    //msg->SetVerboseLevel(1);
    G4AutoDelete::Register(msg);
    fFieldMessenger.Put( msg );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
