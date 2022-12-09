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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 

#include "G4MIRDSmallIntestine.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4EllipticalTube.hh"
#include "G4UnionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDSmallIntestine::Construct(const G4String& volumeName,
						   G4VPhysicalVolume* mother,
						   const G4String& colourName, G4bool wireFrame,G4bool)
{
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
  
  auto* material = new G4HumanPhantomMaterial();
  auto* soft = material -> GetMaterial("soft_tissue");
  delete material;
  
  G4double boxX = 11.*cm;
  G4double boxY = 3.53*cm;
  G4double boxZ = 5*cm;

  auto* smallIntestineBox = new G4Box("smallIntestineBox",boxX,boxY,boxZ);
  
  G4double tubsRmin = 0*cm;
  G4double tubsRmax = 11.*cm;
  G4double tubsZ = 5*cm;
  G4double tubsSphi = 0*degree;
  G4double tubsDphi = 360*degree;
  
  auto* smallIntestineTubs = new G4Tubs("smallIntestineTubs",tubsRmin,tubsRmax,tubsZ,tubsSphi,tubsDphi);

  //G4IntersectionSolid* SmallIntestine = new G4IntersectionSolid("SmallIntestine",smallIntestineTubs,smallIntestineBox,
  auto* filledSmallIntestine1 = new G4IntersectionSolid("filledSmallIntestine1",smallIntestineTubs,smallIntestineBox,
						           nullptr,G4ThreeVector(0*cm,-1.33*cm, 0*cm));

  auto* filledSmallIntestine = new G4IntersectionSolid("filledSmallIntestine",filledSmallIntestine1,smallIntestineTubs,
								      nullptr,G4ThreeVector(0*cm,0.8*cm, 0*cm));

  G4double dx = 2.50*cm; // aU
  G4double dy = 2.50*cm; //bU
  G4double dz = 4.775*cm; //dzU

  auto* AscendingColonUpperLargeIntestine = new G4EllipticalTube("AscendingColon",dx, dy, dz);
 
  dx = 2.50 * cm;//bt
  dy = 1.50 *cm;//ct
  dz = 10.50* cm;//x1t

  auto* TraverseColonUpperLargeIntestine = new G4EllipticalTube("TraverseColon",dx, dy, dz);

  auto* relative_rm =  new G4RotationMatrix();
  relative_rm -> rotateX(90. * degree);
  //relative_rm -> rotateZ(180. * degree);
  relative_rm -> rotateY(90. * degree);
  auto* upperLargeIntestine = new G4UnionSolid("UpperLargeIntestine",
						       AscendingColonUpperLargeIntestine,
						       TraverseColonUpperLargeIntestine,
						       relative_rm, 
						       G4ThreeVector(-8.0 *cm, 0.0*cm,6.275* cm)); //,0,dzU + ct transverse

  dx = 1.88 * cm; //a
  dy = 2.13 *cm; //b
  dz = 7.64 *cm; //(z1-z2)/2
  
  auto* DescendingColonLowerLargeIntestine = new G4EllipticalTube("DiscendingColon",dx, dy, dz);
  
  auto* upperlowerLargeIntestine = new G4UnionSolid("UpperLowerLargeIntestine",
							       upperLargeIntestine,
							       DescendingColonLowerLargeIntestine,
							       nullptr, 
							       G4ThreeVector(-16.72*cm, 0.0*cm,-2.865* cm)); //,0,dzU + ct t
  

  auto* SmallIntestine = new G4SubtractionSolid("SmallIntestine",
							      filledSmallIntestine,
							      upperlowerLargeIntestine,
							      nullptr,
							      G4ThreeVector(8.0*cm,-0.3*cm,-2.775*cm));

  
  auto* logicSmallIntestine = new G4LogicalVolume( SmallIntestine, 
							      soft,
							      "logical"+volumeName,
							      nullptr, nullptr, nullptr);
  auto* rm = new G4RotationMatrix();
  rm->rotateX(180.*degree); 
  rm->rotateY(180.*degree); 
  G4VPhysicalVolume* physSmallIntestine = new G4PVPlacement(rm,     
							    G4ThreeVector(0*cm, -2.66*cm, -13*cm), // Xcheck the spina position the correct placement shuod be this one
							    //G4ThreeVector(0*cm, -5.13*cm, -13*cm), // Xcheck the spina position the correct placement shuod be this one
							    //G4ThreeVector(0*cm, -6*cm, -13*cm),
							    "physical"+volumeName,
							    logicSmallIntestine,
							    mother,
							    false,
							    0, true);
 
  
  // Visualization Attributes
  //G4VisAttributes* SmallIntestineVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* SmallIntestineVisAtt = new G4VisAttributes(colour);
  SmallIntestineVisAtt->SetForceSolid(wireFrame);
  logicSmallIntestine->SetVisAttributes(SmallIntestineVisAtt);
  
  G4cout << "SmallIntestine created !!!!!!" << G4endl;
  
  // Testing SmallIntestine Volume
  G4double SmallIntestineVol = logicSmallIntestine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of SmallIntestine = " << SmallIntestineVol/cm3 << " cm^3" << G4endl;
  
  // Testing SmallIntestine Material
  G4String SmallIntestineMat = logicSmallIntestine->GetMaterial()->GetName();
  G4cout << "Material of SmallIntestine = " << SmallIntestineMat << G4endl;
  
  // Testing Density
  G4double SmallIntestineDensity = logicSmallIntestine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << SmallIntestineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double SmallIntestineMass = (SmallIntestineVol)*SmallIntestineDensity;
  G4cout << "Mass of SmallIntestine = " << SmallIntestineMass/gram << " g" << G4endl;

  return physSmallIntestine;
}
