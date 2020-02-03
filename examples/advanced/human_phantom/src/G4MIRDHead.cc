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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
//
#include "G4MIRDHead.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4EllipticalTube.hh"
#include "G4HumanPhantomColour.hh"

G4MIRDHead::G4MIRDHead()
{

}

G4MIRDHead::~G4MIRDHead()
{

}


G4VPhysicalVolume* G4MIRDHead::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					 const G4String& colourName, G4bool wireFrame, G4bool)
{
  G4cout << "Construct " << volumeName <<" with mother "<<mother->GetName()<<G4endl;
 
  G4HumanPhantomMaterial * material = new G4HumanPhantomMaterial();
  
  G4Material* soft = material -> GetMaterial("soft_tissue");
  
  // MIRD male model
  // Ellipsoid
  G4double ax = 7.0 * cm;
  G4double by = 10.0 * cm;
  G4double cz = 8.50 * cm;
  G4double zcut1 = 0.0 * cm;
  G4double zcut2 = 8.5 * cm;

  G4Ellipsoid* head1 = new G4Ellipsoid("Head1", ax, by, cz, zcut1, zcut2);

  G4double dx = 7.0 * cm;
  G4double dy = 10.0 * cm;
  G4double dz = 7.75 * cm;
 

  G4EllipticalTube* head2 = new G4EllipticalTube("Head2", dx, dy, dz);

  G4UnionSolid* head = new G4UnionSolid("Head",head2,head1,
					0, // Rotation 
					G4ThreeVector(0.* cm, 0.*cm, 7.7500 * cm) );

  G4LogicalVolume* logicHead = new G4LogicalVolume(head, soft,"logical" + volumeName,
						   0, 0,0);
  // Define rotation and position here!
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateX(180.*degree); 
  rm->rotateY(180.*degree); 
  
  G4VPhysicalVolume* physHead = new G4PVPlacement(rm,
						  //G4ThreeVector(0.* cm,0.*cm, 70.75*cm), //FA
						  G4ThreeVector(0.* cm,0.*cm, 77.75*cm),
						  "physicalHead",
						  logicHead,
						  mother,
						  false,
						  0, true);

 
  // Visualization Attributes

  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* HeadVisAtt = new G4VisAttributes(colour);

  HeadVisAtt->SetForceSolid(wireFrame);
  // HeadVisAtt->SetLineWidth(0.7* mm);
  //HeadVisAtt-> SetForceAuxEdgeVisible(true);
  logicHead->SetVisAttributes(HeadVisAtt);

  // Testing Head Volume
  G4double HeadVol = logicHead->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Head = " << HeadVol/cm3 << " cm^3" << G4endl;
  
  // Testing Head Material
  G4String HeadMat = logicHead->GetMaterial()->GetName();
  G4cout << "Material of Head = " << HeadMat << G4endl;
  
  // Testing Density
  G4double HeadDensity = logicHead->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << HeadDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double HeadMass = (HeadVol)*HeadDensity;
  G4cout << "Mass of Head = " << HeadMass/gram << " g" << G4endl;
  
  return physHead;
}
