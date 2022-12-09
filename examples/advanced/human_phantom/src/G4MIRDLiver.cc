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
#include "G4MIRDLiver.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4EllipticalTube.hh"
#include "G4Box.hh"
#include "G4HumanPhantomColour.hh"
#include <cmath>

G4VPhysicalVolume* G4MIRDLiver::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					  const G4String& colourName, G4bool wireFrame, G4bool)
{

  G4cout << "Construct " << volumeName <<" with mother "<<mother->GetName()<<G4endl;
  auto* material = new G4HumanPhantomMaterial();
  auto* soft = material -> GetMaterial("soft_tissue");
  delete material;

  G4double dx= 14.19 *cm; //a
  G4double dy= 7.84 *cm;  //b
  G4double dz= 7.21* cm; //(z2-z1)/2

  auto* firstLiver = new G4EllipticalTube("FirstLiver",dx, dy, dz);

  G4double xx = 20.00 * cm;
  G4double yy = 50.00 * cm;
  G4double zz = 50.00 *cm;

  auto* subtrLiver = new G4Box("SubtrLiver", xx/2., yy/2., zz/2.);

  auto* rm_relative = new G4RotationMatrix();
  rm_relative -> rotateY(32.* degree);
  rm_relative -> rotateZ(40.9* degree);

  //  G4double aa = (1.00/31.51);
  //  G4double bb = (1.00/44.75);
  //  G4double cc = (-1.00/38.76);
  //  G4cout<< aa << " "<< bb << " "<<cc<< G4endl;
 
  //  G4double dd = sqrt(aa*aa + bb*bb + cc*cc);
  //  G4cout<< "dd" << dd << G4endl;
  //  G4cout << aa/dd << "" << bb/dd << " "<< cc/dd <<G4endl;
  //  G4cout << (std::atan(1.42))/deg << G4endl;

  auto* liver = new G4SubtractionSolid("Liver",
 						     firstLiver,subtrLiver,
 						     rm_relative,
 						     G4ThreeVector(10.0*cm,0.0*cm,0.0 *cm));

  auto* logicLiver = new G4LogicalVolume(liver,
 						    soft,
 						    "LiverVolume",
 						    nullptr, nullptr, nullptr);

  // Define rotation and position here!
  auto* rm = new G4RotationMatrix();
  rm->rotateX(180.*degree);
  G4VPhysicalVolume* physLiver = new G4PVPlacement(rm,G4ThreeVector(0. *cm,0. *cm,0.*cm),
						   "physicalLiver",
						   logicLiver,
						   mother,
						   false,
						   0,true);

  // Visualization Attributes
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* LiverVisAtt = new G4VisAttributes(colour);
  LiverVisAtt->SetForceSolid(wireFrame);
  logicLiver->SetVisAttributes(LiverVisAtt);

  G4cout << "Liver created !!!!!!" << G4endl;

  // Testing Liver Volume
  G4double LiverVol = logicLiver->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Liver = " << LiverVol/cm3 << " cm^3" << G4endl;

  // Testing Liver Material
  G4String LiverMat = logicLiver->GetMaterial()->GetName();
  G4cout << "Material of Liver = " << LiverMat << G4endl;

  // Testing Density
  G4double LiverDensity = logicLiver->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LiverDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LiverMass = (LiverVol)*LiverDensity;
  G4cout << "Mass of Liver = " << LiverMass/gram << " g" << G4endl;

  return physLiver;
}
