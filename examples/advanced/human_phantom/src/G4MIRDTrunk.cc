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
#include "G4MIRDTrunk.hh"

#include "G4EllipticalTube.hh"

#include "G4Element.hh"
#include "G4CSGSolid.hh"
#include "G4Material.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"   
#include "G4Colour.hh"


G4MIRDTrunk::G4MIRDTrunk(): physTrunk(0)
{
}

G4MIRDTrunk::~G4MIRDTrunk()
{
}

void G4MIRDTrunk::ConstructTrunk(G4VPhysicalVolume* mother, G4String sex)
{

  // *******************************************************************************//

  // Elements
  G4double A;
  G4double Z;
  A = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);
  A = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",Z = 7.,A);

  // Air Material
  G4double  d = 1.290*mg/cm3;
  G4Material* matAir = new G4Material("Air",d,2);
  matAir->AddElement(elN,0.7);
  matAir->AddElement(elO,0.3);


  // ******************************************************************************//


  G4double Dx = 10.*cm;
  G4double Dy = 20.*cm;
  G4double Dz = 35.*cm;

  G4EllipticalTube* trunk = new G4EllipticalTube("MIRDTrunk", 
						 Dx,
						 Dy,
						 Dz );

  G4LogicalVolume* logicTrunk = new G4LogicalVolume(trunk, 
						    matAir, 
						    "logicalTrunk", 0,0,0);

  physTrunk = new G4PVPlacement(0,G4ThreeVector(),
      				"physicalTrunk",
				logicTrunk,
				mother,
				false,
				0);

  G4VisAttributes* TrunkVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  TrunkVisAtt->SetForceSolid(false);
  logicTrunk->SetVisAttributes(TrunkVisAtt);

 G4cout << "Trunk created !!!!!!" << G4endl;
}
