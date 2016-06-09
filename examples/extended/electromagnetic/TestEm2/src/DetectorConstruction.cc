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
// 
// $Id: DetectorConstruction.cc,v 1.12 2006/10/20 16:03:40 maire Exp $
// GEANT4 tag $Name: geant4-09-01 $

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:nLtot(40),nRtot(50),dLradl(0.5),dRradl(0.1),
 dLlength(0.),dRlength(0.),
 myMaterial(0),magField(0),
 EcalLength(0.),EcalRadius(0.),
 solidEcal(0),logicEcal(0),physiEcal(0)
{
  DefineMaterials();
  SetMaterial("PbWO4");
  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //
  // define few Elements
  //
  G4double a, z;
    
  G4Element* H  = new G4Element("Hydrogen",  "H", z= 1., a=   1.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",  "N", z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,  "O", z= 8., a=  16.00*g/mole);
  G4Element* Ge = new G4Element("Germanium", "Ge",z=32., a=  72.59*g/mole);
  G4Element* W  = new G4Element("Tungsten",  "W", z=74., a= 183.84*g/mole);
  G4Element* Pb = new G4Element("Lead",      "Pb",z=82., a= 207.19*g/mole);
  G4Element* Bi = new G4Element("Bismuth",   "Bi",z=83., a= 208.98*g/mole);

  //
  // define materials
  //
  G4double density;
  G4double fractionmass;  G4int ncomponents, natoms;

  G4Material* Air = 
  new G4Material("Air", density= 1.29*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* H2O = 
  new G4Material("Water", density= 1.00*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);

  new G4Material("Aluminium",   z=13., a= 26.98*g/mole, density= 2.7*g/cm3);

  new G4Material("Iron",        z=26., a= 55.85*g/mole, density= 7.87*g/cm3);
  
  new G4Material("Copper"     , z=29., a= 63.55*g/mole, density= 8.960*g/cm3);
  
  new G4Material("Lead",        z=82., a=207.19*g/mole, density=11.35*g/cm3);
  
  new G4Material("Uranium"    , z=92., a=238.03*g/mole, density= 18.95*g/cm3);
    
  G4Material* BGO = 
  new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
  BGO->AddElement(O , natoms=12);
  BGO->AddElement(Ge, natoms= 3);
  BGO->AddElement(Bi, natoms= 4);

  G4Material* PbWO = 
  new G4Material("PbWO4", density= 8.28*g/cm3, ncomponents=3);
  PbWO->AddElement(O , natoms=4);
  PbWO->AddElement(Pb, natoms=1);
  PbWO->AddElement(W , natoms=1);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4double Radl = myMaterial->GetRadlen();

  dLlength = dLradl*Radl; dRlength = dRradl*Radl;
  EcalLength = nLtot*dLlength;  EcalRadius = nRtot*dRlength;

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //
  // Ecal
  //
  solidEcal = new G4Tubs("Ecal",0.,EcalRadius,0.5*EcalLength,0.,360*deg);
  logicEcal = new G4LogicalVolume( solidEcal,myMaterial,"Ecal",0,0,0);
  physiEcal = new G4PVPlacement(0,G4ThreeVector(),
                                logicEcal,"Ecal",0,false,0);

  G4cout << "Absorber is " << G4BestUnit(EcalLength,"Length")
         << " of " << myMaterial->GetName() << G4endl;
  G4cout << myMaterial << G4endl;     

  //
  //always return the physical World
  //
  return physiEcal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) myMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetLBining(G4ThreeVector Value)
{
  nLtot = (G4int)Value(0);
  if (nLtot > MaxBin) {
    G4cout << "\n ---> warning from SetLBining: "
           << nLtot << " truncated to " << MaxBin << G4endl;
    nLtot = MaxBin;
  }  
  dLradl = Value(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRBining(G4ThreeVector Value)
{
  nRtot = (G4int)Value(0);
  if (nRtot > MaxBin) {
    G4cout << "\n ---> warning from SetRBining: "
           << nRtot << " truncated to " << MaxBin << G4endl;
    nRtot = MaxBin;
  }    
  dRradl = Value(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;		//delete the existing magn field

  if(fieldValue!=0.)			// create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
