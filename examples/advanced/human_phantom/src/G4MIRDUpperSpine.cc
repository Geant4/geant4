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
#include "G4MIRDUpperSpine.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4MIRDUpperSpine::G4MIRDUpperSpine()
{
}

G4MIRDUpperSpine::~G4MIRDUpperSpine()
{
 
}

G4VPhysicalVolume* G4MIRDUpperSpine::ConstructOrgan(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity,
						    G4String volumeName, G4String logicalVolumeName, G4String colourName
						    , G4bool wireFrame )
{
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "Construct " <<volumeName <<" for "<< sex <<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
 
  delete material;

  G4double dx = 2. *cm;
  G4double dy = 2.5 *cm;
  G4double dz = 4.25*cm;

 G4EllipticalTube* upperSpine = new G4EllipticalTube("UpperSpine",dx, dy, dz);

 G4double xx = 20. * cm;
 G4double yy = 10. * cm;
 G4double zz = 5. * cm;

 G4Box* subtraction = new G4Box("box", xx/2., yy/2., zz/2.);

 G4RotationMatrix* matrix = new G4RotationMatrix();
 matrix -> rotateX(-25.* deg); 

 G4SubtractionSolid* upper_spine = new G4SubtractionSolid("upperspine",upperSpine, subtraction,
							  matrix, G4ThreeVector(0., -2.5 * cm, 5.5* cm));

 G4LogicalVolume* logicUpperSpine = new G4LogicalVolume(upper_spine, skeleton, 
							logicalVolumeName,
							0, 0, 0);  
  // Define rotation and position here!
  G4VPhysicalVolume* physUpperSpine = new G4PVPlacement(0,
			        G4ThreeVector(0.0, 5.5 *cm, -3.5 *cm),
      			       "physicalUpperSpine",
  			       logicUpperSpine,
			       mother,
			       false,
			       0, true);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUpperSpine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //G4VisAttributes* UpperSpineVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* UpperSpineVisAtt = new G4VisAttributes(colour);
 
  UpperSpineVisAtt->SetForceSolid(wireFrame);
  logicUpperSpine->SetVisAttributes(UpperSpineVisAtt);

  G4cout << "UpperSpine created !!!!!!" << G4endl;
 
  // Testing UpperSpine Volume
  G4double UpperSpineVol = logicUpperSpine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UpperSpine = " << UpperSpineVol/cm3 << " cm^3" << G4endl;
  
  // Testing UpperSpine Material
  G4String UpperSpineMat = logicUpperSpine->GetMaterial()->GetName();
  G4cout << "Material of UpperSpine = " << UpperSpineMat << G4endl;
  
  // Testing Density
  G4double UpperSpineDensity = logicUpperSpine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UpperSpineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UpperSpineMass = (UpperSpineVol)*UpperSpineDensity;
  G4cout << "Mass of UpperSpine = " << UpperSpineMass/gram << " g" << G4endl;

 
  return physUpperSpine;
}
