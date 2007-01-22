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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4MIRDSpleen.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Ellipsoid.hh"
#include "G4HumanPhantomColour.hh"
G4MIRDSpleen::G4MIRDSpleen()
{
}

G4MIRDSpleen::~G4MIRDSpleen()
{

}

G4VPhysicalVolume* G4MIRDSpleen::ConstructOrgan(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity,G4String volumeName, 
G4String logicalVolumeName, G4String colourName, G4bool wireFrame)
{

  G4cout << "Construct "<< volumeName <<" for " << sex << G4endl;
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax= 3.5 *cm;
 G4double by= 2. *cm;
 G4double cz= 6. * cm; 

 G4Ellipsoid* spleen = new G4Ellipsoid("spleen", ax, by, cz);


  G4LogicalVolume* logicSpleen = new G4LogicalVolume(spleen, soft,
						     logicalVolumeName,
						      0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physSpleen = new G4PVPlacement(0,
						    G4ThreeVector(11. *cm, 3. *cm, 2.*cm), // ztrans = half trunk lenght - z0
      			       "physicalSpleen",
  			       logicSpleen,
			       mother,
			       false,
			       0, true);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicSpleen->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  // G4VisAttributes* SpleenVisAtt = new G4VisAttributes(G4Colour(0.41,0.41,0.41));
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* SpleenVisAtt = new G4VisAttributes(colour);
  SpleenVisAtt->SetForceSolid(wireFrame);
  logicSpleen->SetVisAttributes(SpleenVisAtt);

  G4cout << "Spleen created !!!!!!" << G4endl;

  // Testing Spleen Volume
  G4double SpleenVol = logicSpleen->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Spleen = " << SpleenVol/cm3 << " cm^3" << G4endl;
  
  // Testing Spleen Material
  G4String SpleenMat = logicSpleen->GetMaterial()->GetName();
  G4cout << "Material of Spleen = " << SpleenMat << G4endl;
  
  // Testing Density
  G4double SpleenDensity = logicSpleen->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << SpleenDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double SpleenMass = (SpleenVol)*SpleenDensity;
  G4cout << "Mass of Spleen = " << SpleenMass/gram << " g" << G4endl;


  
  return physSpleen;
}
