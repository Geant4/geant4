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
#include "G4PhantomBuilder.hh"
#include "G4VBodyFactory.hh"
#include "G4MIRDBodyFactory.hh"
#include "G4ORNLFemaleBodyFactory.hh"
#include "G4ORNLMaleBodyFactory.hh"
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4PhantomBuilder::G4PhantomBuilder(): model("MIRD")
{  
  // sex can be "female" or "male"
  body = 0;
  motherVolume = 0;
  headVolume = 0;
  trunkVolume = 0;
  leftLegVolume =0;
  rightLegVolume =0;
  maleGenitaliaVolume = 0;  
}

G4PhantomBuilder::~G4PhantomBuilder()
{
} 
void G4PhantomBuilder::BuildTrunk(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("G4PhantomBuilder::BuildTrunk()", "human_phantom0014", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  trunkVolume = body -> CreateOrgan("Trunk", motherVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftLeg(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftLeg()", "human_phantom0015", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  leftLegVolume = body -> CreateOrgan("LeftLeg", motherVolume, colourName, solidVis, sensitivity);
}
void G4PhantomBuilder::BuildRightLeg(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRightLeg()", "human_phantom0016", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  rightLegVolume = body -> CreateOrgan("RightLeg", motherVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftLegBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (leftLegVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftLegBone()", "human_phantom0017", FatalException, "The left leg volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  leftLegVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("LeftLegBone", leftLegVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightLegBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRightLegBone()", "human_phantom0018", FatalException, "The right leg volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << rightLegVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("RightLegBone", rightLegVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftArmBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftArmBone()", "human_phantom0019", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("LeftArmBone" ,trunkVolume,colourName,solidVis, sensitivity);
}
void G4PhantomBuilder::BuildRightArmBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRightArmBone()", "human_phantom0020", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("RightArmBone",trunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftScapula(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftScapula()", "human_phantom0021", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("LeftScapula",trunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightScapula(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRightScapula()", "human_phantom0022", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("RightScapula",trunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftClavicle(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftClavicle()", "human_phantom0023", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("LeftClavicle",trunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightClavicle(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRightClavicle()", "human_phantom0024", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("RightClavicle",trunkVolume,colourName,solidVis, sensitivity);
}


void G4PhantomBuilder::BuildHead(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("G4PhantomBuilder::BuildHead()", "human_phantom0025", FatalException, "The mother volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  headVolume = body -> CreateOrgan("Head",motherVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildSkull(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (headVolume == 0)
    G4Exception("G4PhantomBuilder::BuildSkull()", "human_phantom0026", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  headVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan( "Skull",headVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildUpperSpine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (headVolume == 0)
    G4Exception("G4PhantomBuilder::BuildUpperSpine()", "human_phantom0027", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  headVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("UpperSpine",headVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildMiddleLowerSpine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildMiddleLowerSpine()", "human_phantom0028", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("MiddleLowerSpine",trunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildPelvis(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildPelvis()", "human_phantom0029", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan( "Pelvis",trunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildBrain(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (headVolume == 0)
    G4Exception("G4PhantomBuilder::BuildBrain()", "human_phantom0030", FatalException, "The head volume is missing !!!!!");

  body -> CreateOrgan("Brain",headVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildHeart(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildHeart()", "human_phantom0031", FatalException, "The trunk volume is missing !!!!!");
  body -> CreateOrgan("Heart", trunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftLung(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftLung()", "human_phantom0032", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("LeftLung",trunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightLung(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRightLung()", "human_phantom0033", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("RightLung",trunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildStomach(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildStomach()", "human_phantom0034", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("Stomach",trunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRibCage(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRibCage()", "human_phantom0035", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("RibCage",trunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildSpleen(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildSpleen()", "human_phantom0036", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("Spleen", trunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildUpperLargeIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildUpperLargeIntestine()", "human_phantom0037", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("UpperLargeIntestine",trunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLowerLargeIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLowerLargeIntestine()", "human_phantom0038", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("LowerLargeIntestine", trunkVolume, colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildSmallIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildSamllIntestine()", "human_phantom0039", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("SmallIntestine",trunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftKidney(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftKidney()", "human_phantom0040", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("LeftKidney", trunkVolume,colourName, solidVis, sensitivity);
}
void G4PhantomBuilder::BuildRightKidney(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRightKidney()", "human_phantom0041", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("RightKidney",trunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftAdrenal(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftAdrenal()", "human_phantom0042", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("LeftAdrenal", trunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightAdrenal(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildRightAdrenal()", "human_phantom0043", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("RightAdrenal", trunkVolume,colourName, solidVis, sensitivity);
}


void G4PhantomBuilder::BuildLiver(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLiver()", "human_phantom0044", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("Liver", trunkVolume,colourName, solidVis, sensitivity);
}
void G4PhantomBuilder::BuildPancreas(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildPancreas()", "human_phantom0045", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("Pancreas",trunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildUrinaryBladder(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildUrinaryBladder()", "human_phantom0046", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("UrinaryBladder",trunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildThyroid(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (headVolume == 0)
    G4Exception("G4PhantomBuilder::BuildThyroid()", "human_phantom0047", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("Thyroid",headVolume, colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildThymus(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildThymus()", "human_phantom0048", FatalException, "The trunk volume is missing !!!!!");
  
  body -> CreateOrgan("Thymus",trunkVolume, colourName,solidVis, sensitivity);
}


G4VPhysicalVolume* G4PhantomBuilder::GetPhantom()
{
  return motherVolume;
}

void G4PhantomBuilder::SetMotherVolume(G4VPhysicalVolume* mother)
{
  motherVolume = mother;
}


void G4PhantomBuilder::SetModel(G4String modelFlag)
{
  model = modelFlag;

  if(model=="MIRD" || model =="MIX") body = new G4MIRDBodyFactory();
  if(model=="ORNLFemale")
  {

#ifdef G4LIB_USE_GDML
    body = new G4ORNLFemaleBodyFactory();
#else
    G4cout << model << " Working with GDML only! set G4LIB_USE_GDML 1" << G4endl;
#endif
  }

  if(model=="ORNLMale") 
  {
#ifdef G4LIB_USE_GDML
    body = new G4ORNLMaleBodyFactory();
#else
    G4cout << model << " Working with GDML only! set G4LIB_USE_GDML 1" << G4endl;
#endif
  }
}
