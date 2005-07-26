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
#include "G4PhantomBuilder.hh"
#include "G4VBodyFactory.hh"
#include "G4MIRDBodyFactory.hh"
#include "G4ORNLBodyFactory.hh"
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

G4PhantomBuilder::G4PhantomBuilder(): sex("female")
{  
  // sex can be "female" or "male"
  motherVolume = 0;
  headVolume = 0;
  trunkVolume = 0;
  
  body = new G4MIRDBodyFactory();
}

G4PhantomBuilder::~G4PhantomBuilder()
{
  delete body;
} 

void G4PhantomBuilder::BuildWorld()
{
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

  // World Volume
  G4double worldSize = 1.*m ;

  G4Box* world = new G4Box("world", worldSize, worldSize, worldSize);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(world, 
						    matAir, 
						    "logicalWorld", 0, 0,0);

  motherVolume = new G4PVPlacement(0,G4ThreeVector(),
						"physicalWorld",
						logicWorld,
						0,
						false,
						0);
}
void G4PhantomBuilder::BuildHead()
{ 
  if (motherVolume == 0)
    G4Exception("The world volume is missing !!!!!");

  headVolume = body -> CreateHead(motherVolume,sex);
}

void G4PhantomBuilder::BuildTrunk()
{ 
  if (motherVolume == 0)
    G4Exception("The world volume is missing !!!!!");

  trunkVolume = body -> CreateTrunk(motherVolume,sex);
}
//void G4PhantomBuilder::BuildBrain()
//{ 
  // if (headVolume == 0)
  //  G4Exception("The head volume is missing !!!!!");

  // body -> CreateBrain(headVolume);
//}

G4VPhysicalVolume* G4PhantomBuilder::GetPhantom()
{
  return motherVolume;
}
