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
#include "G4Colour.hh"
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
void G4PhantomBuilder::BuildTrunk(G4bool sensitivity, G4String volumeName, G4String logicalVolumeName, G4String colourName, G4bool solidVis)
{ 
  if (motherVolume == 0)
    G4Exception("The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  trunkVolume = body -> CreateOrgan(motherVolume,sex, sensitivity, volumeName, logicalVolumeName, colourName, solidVis);
}

void G4PhantomBuilder::BuildHead(G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  headVolume = body -> CreateOrgan(motherVolume,sex, sensitivity, "Head", "Default", "pink", false);
}

void G4PhantomBuilder::BuildPelvis(G4bool sensitivity)
{ 
  if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity, "Pelvis","Default", "pink", false);
}
/*
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
*/

void G4PhantomBuilder::BuildBrain(G4bool sensitivity)
{ 
 if (headVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(headVolume,sex,sensitivity, "Brain", "Default", "pink", false);
}

void G4PhantomBuilder::BuildHeart(G4bool sensitivity)
{ 
  //  if (trunkVolume == 0)
//    G4Exception("The trunk volume is missing !!!!!");

//    body -> CreateHeart(trunkVolume,sex,sensitivity);
}

void G4PhantomBuilder::BuildLeftLung(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
    G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity, "LeftLung", "Default", "pink", false);
}

void G4PhantomBuilder::BuildRightLung(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
    G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity, "RightLung", "Default", "pink", false);
}

void G4PhantomBuilder::BuildStomach(G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity, "Stomach", "Default", "pink", false);
}

void G4PhantomBuilder::BuildUpperLargeIntestine(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity, "UpperLargeIntestine", "Default", "pink", false);
}

void G4PhantomBuilder::BuildLowerLargeIntestine(G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("The trunk volume is missing !!!!!");

   body -> CreateOrgan(trunkVolume,sex,sensitivity, "LowerLargeIntestine", "Default", "pink",false );
}
/*
void G4PhantomBuilder::BuildEsophagus(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateEsophagus(trunkVolume,sex,sensitivity);
}
*/
void G4PhantomBuilder::BuildLeftKidney(G4bool sensitivity)
{ 
 if (trunkVolume == 0)
    G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity,"LeftKidney", "Default", "pink", false);
}
void G4PhantomBuilder::BuildRightKidney(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
    G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity, "RightKidney", "Default", "pink", false);
}

/*
void G4PhantomBuilder::BuildAdrenal(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

   body -> CreateAdrenal(trunkVolume,sex,sensitivity);
}
*/
void G4PhantomBuilder::BuildLiver(G4bool sensitivity)
{ 
   // if (trunkVolume == 0)
//    G4Exception("The trunk volume is missing !!!!!");

//    body -> CreateLiver(trunkVolume,sex,sensitivity);
}
void G4PhantomBuilder::BuildPancreas(G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity,"Pancreas", "Default", "pink", false);
}

void G4PhantomBuilder::BuildUrinaryBladder(G4bool sensitivity)
{ 
  if (trunkVolume == 0)
   G4Exception("The trunk volume is missing !!!!!");

    body -> CreateOrgan(trunkVolume,sex,sensitivity, "UrinaryBladder", "Default", "pink", false);
}

void G4PhantomBuilder::BuildThyroid(G4bool sensitivity)
{ 
   // if (headVolume == 0)
//    G4Exception("The trunk volume is missing !!!!!");

//    body -> CreateThyroid(headVolume,sex,sensitivity);
}


G4VPhysicalVolume* G4PhantomBuilder::GetPhantom()
{
  return motherVolume;
}

void G4PhantomBuilder::SetMotherVolume(G4VPhysicalVolume* mother)
{
  motherVolume = mother;
}

void G4PhantomBuilder::SetSex(G4String sexFlag)
{
  sex = sexFlag;
}

void G4PhantomBuilder::SetModel(G4String modelFlag)
{
  model = modelFlag;

  if(model=="MIRD") body = new G4MIRDBodyFactory();
  if(model=="ORNL") body = new G4ORNLBodyFactory();
}
