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

#include "G4VisAttributes.hh"

G4PhantomBuilder::G4PhantomBuilder(): sex("Female"), model("MIRD")
{  
  // sex can be "female" or "male"
  motherVolume = 0;
  headVolume = 0;
  trunkVolume = 0;
  neckVolume = 0;
  maleGenitaliaVolume = 0;  
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

  // Visualization Attributes
  G4VisAttributes* WorldVisAtt = new G4VisAttributes(G4VisAttributes::Invisible);
  WorldVisAtt->SetForceSolid(false);
  logicWorld->SetVisAttributes(WorldVisAtt);


}

void G4PhantomBuilder::BuildHead(G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("The world volume is missing !!!!!");

  headVolume = body -> CreateHead(motherVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildTrunk(G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("The world volume is missing !!!!!");

  trunkVolume = body -> CreateTrunk(motherVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildLegs(G4bool sensitivity)
{ 
   if (motherVolume == 0)
   G4Exception("The world volume is missing !!!!!");

   legsVolume = body -> CreateLegs(motherVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildNeck(G4bool sensitivity)
{ 
   if (motherVolume == 0)
   G4Exception("The world volume is missing !!!!!");

   neckVolume = body -> CreateNeck(motherVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildUpperSpine(G4bool sensitivity)
{ 
   if (headVolume == 0)
   G4Exception("The head volume is missing !!!!!");

   body -> CreateUpperSpine(headVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildMiddleLowerSpine(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateMiddleLowerSpine(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildSpleen(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateSpleen(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildLegBone(G4bool sensitivity)
{ 
   if (legsVolume == 0)
   G4Exception("The legs volume volume is missing !!!!!");

   body -> CreateLegBone(legsVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildArmBone(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateArmBone(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildSkull(G4bool sensitivity)
{ 
   if (headVolume == 0)
   G4Exception("The head volume is missing !!!!!");

   body -> CreateSkull(headVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildRibCage(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateRibCage(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildPelvis(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreatePelvis(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildScapulae(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateScapulae(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildClavicles(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateClavicles(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildBreast(G4bool sensitivity)
{ 
   if (motherVolume == 0)
   G4Exception("The world volume is missing !!!!!");

   body -> CreateBreast(motherVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildUterus(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateUterus(trunkVolume,sex,sensitivity);
}
void G4PhantomBuilder::BuildOvary(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateOvary(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildMaleGenitalia(G4bool sensitivity)
{ 
   if (motherVolume == 0)
   G4Exception("The world volume is missing !!!!!");

   maleGenitaliaVolume = body -> CreateMaleGenitalia(motherVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildTestes(G4bool sensitivity)
{ 
   if (maleGenitaliaVolume == 0)
   G4Exception("The world volume is missing !!!!!");

   body -> CreateTestes(maleGenitaliaVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildBrain(G4bool sensitivity)
{ 
   if (headVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateBrain(headVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildHeart(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateHeart(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildLung(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateLung(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildStomach(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateStomach(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildUpperLargeIntestine(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateUpperLargeIntestine(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildLowerLargeIntestine(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateLowerLargeIntestine(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildEsophagus(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateEsophagus(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildKidney(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateKidney(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildAdrenal(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateAdrenal(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildLiver(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateLiver(trunkVolume,sex,sensitivity);
}
void G4PhantomBuilder::BuildPancreas(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreatePancreas(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildUrinaryBladder(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateUrinaryBladder(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildThyroid(G4bool sensitivity)
{ 
   if (headVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateThyroid(headVolume,sex,sensitivity);
}


G4VPhysicalVolume* G4PhantomBuilder::GetPhantom()
{
  return motherVolume;
}

void G4PhantomBuilder::SetSex(G4String sexFlag)
{
  sex = sexFlag;
}

void G4PhantomBuilder::SetModel(G4String modelFlag)
{
  model = modelFlag;

  if(model=="MIRD")
    body = new G4MIRDBodyFactory();
  if(model=="ORNL")
    body = new G4ORNLBodyFactory();
}
