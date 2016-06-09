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
/// \file electromagnetic/TestEm11/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fDefaultMaterial(0),fPhysiWorld(0),fMagField(0)
{
  // default parameter values of the absorbers
  fNbOfAbsor = 1;
  fAbsorThickness[0] = 0*mm;        //dummy, for initialization   
  fAbsorThickness[1] = 1*mm;  
  fAbsorSizeYZ       = 1.*mm;
  for (G4int iAbs=0; iAbs<MaxAbsor; iAbs++) {
    fNbOfDivisions[iAbs]  = 1;
  }  
  ComputeParameters();

  // materials
  DefineMaterials();
  SetAbsorMaterial(1,"G4_Si");

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  
  G4bool isotopes = false;
  
  man->FindOrBuildMaterial("G4_Al", isotopes);
  man->FindOrBuildMaterial("G4_Si", isotopes);
  man->FindOrBuildMaterial("G4_Fe", isotopes);
  man->FindOrBuildMaterial("G4_Cu", isotopes);  
  man->FindOrBuildMaterial("G4_Ge", isotopes);
  man->FindOrBuildMaterial("G4_Mo", isotopes);
  man->FindOrBuildMaterial("G4_Ta", isotopes);
  man->FindOrBuildMaterial("G4_W" , isotopes);
  man->FindOrBuildMaterial("G4_Au", isotopes);
  man->FindOrBuildMaterial("G4_Pb", isotopes);  
  man->FindOrBuildMaterial("G4_PbWO4", isotopes);
  man->FindOrBuildMaterial("G4_SODIUM_IODIDE", isotopes);
  
  man->FindOrBuildMaterial("G4_AIR"  , isotopes);
  man->FindOrBuildMaterial("G4_WATER", isotopes);
  
  G4Element* H = man->FindOrBuildElement("H", isotopes); 
  G4Element* O = man->FindOrBuildElement("O", isotopes);
  
  G4Material* H2O = 
  new G4Material("Water", 1.000*g/cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  G4double density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4Material* Galactic =   
  new G4Material("Galactic", 1., 1.008*g/mole, density,
                             kStateGas,temperature,pressure);
                             
  fDefaultMaterial = Galactic;
  
  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeParameters()
{
  // Compute total thickness of absorbers
  fAbsorSizeX = 0.;
  for (G4int iAbs=1; iAbs<=fNbOfAbsor; iAbs++) {
    fAbsorSizeX += fAbsorThickness[iAbs];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // complete the Calor parameters definition
  ComputeParameters();

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //
  // World
  //
  G4Box* solidWorld =
    new G4Box("World",                                                //name
               fAbsorSizeX/2,fAbsorSizeYZ/2,fAbsorSizeYZ/2);        //size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,                //solid
                        fDefaultMaterial,        //material
                        "World");                //name

  fPhysiWorld = 
    new G4PVPlacement(0,                        //no rotation
                        G4ThreeVector(),                //at (0,0,0)
                      logicWorld,                //logical volume
                      "World",                        //name
                       0,                        //mother volume
                       false,                        //no boolean operation
                       0);                       //copy number
                                 
  //
  // Absorbers
  //
  fXfront[0] = -0.5*fAbsorSizeX;
  //
  for (G4int k=1; k<=fNbOfAbsor; k++) {
    G4Material* material = fAbsorMaterial[k];
    G4String matname = material->GetName();
      
    G4Box* solidAbsor =
      new G4Box(matname,fAbsorThickness[k]/2,fAbsorSizeYZ/2,fAbsorSizeYZ/2);

    G4LogicalVolume* logicAbsor =
      new G4LogicalVolume(solidAbsor,            // solid
                          material,                 // material
                          matname);             // name
                                     
    fXfront[k] = fXfront[k-1] + fAbsorThickness[k-1];    
    G4double xcenter = fXfront[k]+0.5*fAbsorThickness[k];
    G4ThreeVector position = G4ThreeVector(xcenter,0.,0.);
  
      new G4PVPlacement(0,                           //no rotation
                              position,                   //position
                        logicAbsor,             //logical volume        
                        matname,                   //name
                        logicWorld,                //mother
                        false,                     //no boulean operat
                        k);                       //copy number
    
    // divisions, if any
    //
    G4double LayerThickness = fAbsorThickness[k]/fNbOfDivisions[k];
    G4Box* solidLayer =   
      new G4Box(matname,LayerThickness/2,fAbsorSizeYZ/2,fAbsorSizeYZ/2);
                       
    G4LogicalVolume* logicLayer =
      new G4LogicalVolume(solidLayer,        //solid
                          material,        //material
                          matname);         //name
                                       
      new G4PVReplica(matname,                        //name
                            logicLayer,                //logical volume
                            logicAbsor,                //mother
                      kXAxis,                        //axis of replication
                      fNbOfDivisions[k],                //number of replica
                      LayerThickness);                //witdth of replica    
  }


  PrintParameters();

  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The Absorber is " << fNbOfAbsor << " layers of:";
  for (G4int i=1; i<=fNbOfAbsor; i++)
     {
      G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() <<": "
              << std::setw(6) << G4BestUnit(fAbsorThickness[i],"Length")
              << "  divided in " << fNbOfDivisions[i] << " slices";
     }
  G4cout << "\n-------------------------------------------------------------\n"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (MaxAbsor-1))
    { G4cout << "\n ---> warning from SetfNbOfAbsor: "
             << ival << " must be at least 1 and and most " << MaxAbsor-1
             << ". Command refused" << G4endl;
      return;
    }
  fNbOfAbsor = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4int iabs, const G4String& material)
{
  // search the material by its name
  //
  if (iabs > fNbOfAbsor || iabs <= 0)
    { G4cout << "\n --->warning from SetfAbsorMaterial: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }

  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) fAbsorMaterial[iabs] = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorThickness(G4int iabs,G4double val)
{
  // change Absorber thickness
  //
  if (iabs > fNbOfAbsor || iabs <= 0)
    { G4cout << "\n --->warning from SetfAbsorThickness: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorThickness: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorThickness[iabs] = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorSizeYZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorSizeYZ: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorSizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfDivisions(G4int iabs, G4int ival)
{
  // set the number of divisions
  //
  if (iabs > fNbOfAbsor || iabs < 1)
    { G4cout << "\n --->warning from SetNbOfDivisions: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }
      
  if (ival < 1)
    { G4cout << "\n --->warning from SetNbOfDivisions: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fNbOfDivisions[iabs] = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  //
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(fMagField) delete fMagField;                //delete the existing magn field

  if(fieldValue!=0.)                        // create a new one if non nul
  { fMagField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(fMagField);
    fieldMgr->CreateChordFinder(fMagField);
  } else {
    fMagField = 0;
    fieldMgr->SetDetectorField(fMagField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
