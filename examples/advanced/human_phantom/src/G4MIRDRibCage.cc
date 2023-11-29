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
//
#include "G4MIRDRibCage.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4SubtractionSolid.hh"
#include "G4EllipticalTube.hh"
#include "G4PVReplica.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4HumanPhantomColour.hh"
G4MIRDRibCage::G4MIRDRibCage():fPhysRib1(nullptr), fPhysRib2(nullptr), fPhysRib3(nullptr),
					    fPhysRib4(nullptr), fPhysRib5(nullptr), fPhysRib6(nullptr),
					    fPhysRib7(nullptr), fPhysRib8(nullptr), fPhysRib9(nullptr),
					    fPhysRib10(nullptr), fPhysRib11(nullptr), fPhysRib12(nullptr)
{;}

G4VPhysicalVolume* G4MIRDRibCage::Construct(const G4String& volumeName, G4VPhysicalVolume* mother,
					    const G4String& colourName, G4bool wireFrame, G4bool)
					    
{
  auto* material = new G4HumanPhantomMaterial();
   
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
   
  auto* skeleton = material -> GetMaterial("skeleton");
  auto* soft = material -> GetMaterial("soft_tissue");
 
  delete material; 

  G4double dx= 17. *cm; // a2
  G4double dy= 9.8 * cm; //b2
  G4double thickness= 32.4 * cm; // z2/2 of cage

  auto* outCage = new G4EllipticalTube("outCage",dx, dy, thickness/2.);

  dx = 16.4 * cm; // a1
  dy = 9.2 * cm; // b1
  G4double dz = 34. *cm; // z2/2

  auto* inCage = new G4EllipticalTube("inCage",dx, dy, dz/2.);

  auto* cage = new G4SubtractionSolid("Cage",
						    outCage,
						    inCage, nullptr, G4ThreeVector(0.*cm, 0.*cm, 0. * cm)); 

 
  auto* logicRibCage = new G4LogicalVolume(cage, soft, "logicalCage", nullptr, nullptr, nullptr);

  G4VPhysicalVolume* physRibCage = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, thickness/2. + 0.1 * cm),
						     // with respect to the trunk
						     "physicalRibCage",
						     logicRibCage,
						     mother,
						     false,
						     0, true);
	
  G4double xx = 17.*cm;
  G4double yy = 9.8*cm;
  G4double ribThickness = 1.4*cm;
  auto* rib_out = new G4EllipticalTube("rib_out",xx, yy, ribThickness/2.);	
  
  xx = 16.5 *cm;
  yy = 9.3 * cm;
  G4double zz = 1.5 * cm;  
  auto* rib_in = new G4EllipticalTube("rib_in",xx, yy, zz/2.);
  auto* rib = new G4SubtractionSolid("rib",rib_out, rib_in);

  auto* logicRib= new G4LogicalVolume(rib, skeleton, "logical" + volumeName, nullptr, nullptr, nullptr);

  fPhysRib1 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (- 32.2*cm/2. + 0.8 *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib2 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, ( - 32.2*cm/2. + 0.8 *cm + 2.8 *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib3 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 * cm + 5.6 *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib4 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 * cm + 8.4 *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib5 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 * cm + 11.2 *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib6 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. +  0.8 * cm + 14. *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib7 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 *cm + 16.8 *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib8 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 *cm + 19.6 *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib9 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8*cm + 22.4 *cm)),
			       // with respect to the trunk
			       "physicalRib",
			       logicRib,
			       physRibCage,
			       false,
			       0, true);

  fPhysRib10 = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8*cm + 25.2 *cm)),
				// with respect to the trunk
				"physicalRib",
				logicRib,
				physRibCage,
				false,
				0, true);

  fPhysRib11 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8*cm + 28. *cm)),
				// with respect to the trunk
				"physicalRib",
				logicRib,
				physRibCage,
				false,
				0, true);

  fPhysRib12 = new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8*cm + 30.8 *cm)),
				// with respect to the trunk
				"physicalRib",
				logicRib,
				physRibCage,
				false,
				0, true);
  
  // Visualization Attributes
  logicRibCage -> SetVisAttributes(G4VisAttributes::GetInvisible());

  //G4VisAttributes* RibCageVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));

  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* RibCageVisAtt = new G4VisAttributes(colour); 
  RibCageVisAtt->SetForceSolid(wireFrame);
  logicRib->SetVisAttributes(RibCageVisAtt);

  G4cout << "RibCage created !!!!!!" << G4endl;
  // Testing Pelvis Volume
  G4double RibCageVol = logicRib->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RibCage = " << ((RibCageVol)*12.)/cm3 << " cm^3" << G4endl;
  
  // Testing RibCage Material
  G4String RibCageMat = logicRib->GetMaterial()->GetName();
  G4cout << "Material of RibCage = " << RibCageMat << G4endl;
  
  // Testing Density
  G4double RibCageDensity = logicRib->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RibCageDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RibCageMass = (RibCageVol)* RibCageDensity * 12;// 12 is the total number of ribs;
  G4cout << "Mass of RibCage = " << (RibCageMass)/gram << " g" << G4endl;

  return physRibCage;
}
