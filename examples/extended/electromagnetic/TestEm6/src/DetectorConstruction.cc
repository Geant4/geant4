//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

//
// $Id: DetectorConstruction.cc,v 1.3 2003-01-31 12:54:23 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:pBox(0), lBox(0), BoxSize(500*m),
 aMaterial(0), 
 magField (0)
{
  // create UserLimits
  userLimits = new G4UserLimits();

  // create commands for interactive definition of the detector  
  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
 G4double a, z, density; 

 z = 4.;
 a = 9.012182*g/mole;
 density = 1.848*g/cm3;
 G4Material* Be = new G4Material("Beryllium", z, a, density);
 
 z = 6.;
 a = 12.011*g/mole;
 density = 2.265*g/cm3;
 G4Material* C = new G4Material("Carbon", z, a, density);

 z = 26.;
 a = 55.85*g/mole;
 density = 7.870*g/cm3;
 G4Material* Fe  = new G4Material("Iron"    , z, a, density);

 G4cout << *(G4Material::GetMaterialTable()) << G4endl;

 //default materials of the calorimeter
 aMaterial = Fe;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4Box*
  sBox = new G4Box("Container",				//its name
                   BoxSize/2,BoxSize/2,BoxSize/2);	//its dimensions
		   			                      
  lBox = new G4LogicalVolume(sBox,			//its shape
                             aMaterial,			//its material
                             "Container");		//its name
			     
  lBox->SetUserLimits(userLimits);			     
                                   
  pBox = new G4PVPlacement(0,				//no rotation
  			   G4ThreeVector(),		//at (0,0,0)
                           "Container",			//its name
                           lBox,			//its logical volume
                           0,	       		        //its mother  volume
                           false,			//no boolean operation
                           0);				//copy number
                           
  //                                        
  // Visualization attributes
  //

  G4VisAttributes* visAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAtt->SetVisibility(true);
  lBox->SetVisAttributes(visAtt);
  
  //
  //always return the root volume
  //
  
  PrintParameters();
  return pBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(BoxSize,"Length")
         << " of " << aMaterial->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
     {aMaterial = pttoMaterial;
      lBox->SetMaterial(aMaterial); 
     }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  BoxSize = value;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr 
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if (magField) delete magField;	//delete the existing magn field
  
  if (fieldValue!=0.)			// create a new one if non nul
    {
      magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));        
      fieldMgr->SetDetectorField(magField);
      fieldMgr->CreateChordFinder(magField);
    }
   else
    {
      magField = 0;
      fieldMgr->SetDetectorField(magField);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStepSize(G4double val)
{
  // set the maximum allowed step size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetMaxStepSize: maxStep " 
             << val  << " out of range. Command refused" << G4endl;
      return;
    }       
  userLimits->SetMaxAllowedStep(val);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
