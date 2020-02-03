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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
//
#include "G4MIRDLeftLeg.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4Cons.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4MIRDLeftLeg::G4MIRDLeftLeg()
{
}

G4MIRDLeftLeg::~G4MIRDLeftLeg()
{
}


G4VPhysicalVolume* G4MIRDLeftLeg::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
					    const G4String& colourName, G4bool wireFrame,G4bool)
{
 
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
  G4Material* soft = material -> GetMaterial("soft_tissue");
 
  G4double rmin1 = 0.* cm;
  G4double rmin2 = 0.* cm;
  G4double dz= 80.0 * cm; 
  G4double rmax1= 2.0 * cm;
  G4double rmax2= 10. * cm;
  G4double startphi= 0.* degree;
  G4double deltaphi= 360. * degree;

  G4Cons* leg1 = new G4Cons("Leg1",  
			    rmin1, rmax1, 
			    rmin2, rmax2, dz/2., 
			    startphi, deltaphi);
  
  G4LogicalVolume* logicLeftLeg = new G4LogicalVolume(leg1,
						      soft,
						      "logical" + volumeName,
						      0, 0, 0);
						   
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateX(180.*degree); 
  rm->rotateY(180.*degree);
  G4VPhysicalVolume* physLeftLeg = new G4PVPlacement(rm,
						     //G4ThreeVector(10. * cm, 0. * cm, -47. *cm), //FA
						     G4ThreeVector(10. * cm, 0. * cm, -40. *cm),
						     "physicalLeftLeg",
						     logicLeftLeg,
						     mother,
						     false,
						     0, true);

 
  // Visualization Attributes
  //G4VisAttributes* LeftLegVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* LeftLegVisAtt = new G4VisAttributes(colour);
  LeftLegVisAtt->SetForceSolid(wireFrame);
  logicLeftLeg->SetVisAttributes(LeftLegVisAtt);

  G4cout << "LeftLeg created !!!!!!" << G4endl;

  // Testing LeftLeg Volume
  G4double LeftLegVol = logicLeftLeg->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LeftLeg = " << LeftLegVol/cm3 << " cm^3" << G4endl;
  
  // Testing LeftLeg Material
  G4String LeftLegMat = logicLeftLeg->GetMaterial()->GetName();
  G4cout << "Material of LeftLeg = " << LeftLegMat << G4endl;
  
  // Testing Density
  G4double LeftLegDensity = logicLeftLeg->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftLegDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftLegMass = (LeftLegVol)*LeftLegDensity;
  G4cout << "Mass of LeftLeg = " << LeftLegMass/gram << " g" << G4endl;

  
  return physLeftLeg;
}
