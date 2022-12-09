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

G4PhantomBuilder::G4PhantomBuilder(): fModel("MIRD")
{  
  // sex can be "female" or "male"
  fBody = nullptr;
  fMotherVolume = nullptr;
  fHeadVolume = nullptr;
  fTrunkVolume = nullptr;
  fLeftLegVolume =nullptr;
  fRightLegVolume = nullptr;
  fMaleGenitaliaVolume = nullptr;  
}

void G4PhantomBuilder::BuildTrunk(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fMotherVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildTrunk()", "human_phantom0014", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  fMotherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fTrunkVolume = fBody -> CreateOrgan("Trunk", fMotherVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftLeg(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fMotherVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLeftLeg()", "human_phantom0015", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  fMotherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fLeftLegVolume = fBody -> CreateOrgan("LeftLeg", fMotherVolume, colourName, solidVis, sensitivity);
}
void G4PhantomBuilder::BuildRightLeg(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fMotherVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRightLeg()", "human_phantom0016", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  fMotherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fRightLegVolume = fBody -> CreateOrgan("RightLeg", fMotherVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftLegBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fLeftLegVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLeftLegBone()", "human_phantom0017", FatalException, "The left leg volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  fLeftLegVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("LeftLegBone", fLeftLegVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightLegBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRightLegBone()", "human_phantom0018", FatalException, "The right leg volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fRightLegVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("RightLegBone", fRightLegVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftArmBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLeftArmBone()", "human_phantom0019", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fTrunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("LeftArmBone", fTrunkVolume, colourName, solidVis, sensitivity);
}
void G4PhantomBuilder::BuildRightArmBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRightArmBone()", "human_phantom0020", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fTrunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("RightArmBone", fTrunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftScapula(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLeftScapula()", "human_phantom0021", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fTrunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("LeftScapula",fTrunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightScapula(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRightScapula()", "human_phantom0022", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  fTrunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("RightScapula",fTrunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftClavicle(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLeftClavicle()", "human_phantom0023", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  fTrunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("LeftClavicle", fTrunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightClavicle(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRightClavicle()", "human_phantom0024", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fTrunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("RightClavicle",fTrunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildHead(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fMotherVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildHead()", "human_phantom0025", FatalException, "The mother volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  fMotherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fHeadVolume = fBody -> CreateOrgan("Head", fMotherVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildSkull(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fHeadVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildSkull()", "human_phantom0026", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fHeadVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan( "Skull", fHeadVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildUpperSpine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fHeadVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildUpperSpine()", "human_phantom0027", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fHeadVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("UpperSpine", fHeadVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildMiddleLowerSpine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildMiddleLowerSpine()", "human_phantom0028", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fTrunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan("MiddleLowerSpine",fTrunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildPelvis(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildPelvis()", "human_phantom0029", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan( "Pelvis", fTrunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildBrain(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fHeadVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildBrain()", "human_phantom0030", FatalException, "The head volume is missing !!!!!");

  fBody -> CreateOrgan("Brain", fHeadVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildHeart(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildHeart()", "human_phantom0031", FatalException, "The trunk volume is missing !!!!!");
  fBody -> CreateOrgan("Heart", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftLung(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLeftLung()", "human_phantom0032", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("LeftLung", fTrunkVolume,colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightLung(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRightLung()", "human_phantom0033", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("RightLung", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildStomach(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildStomach()", "human_phantom0034", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("Stomach", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRibCage(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRibCage()", "human_phantom0035", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("RibCage", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildSpleen(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildSpleen()", "human_phantom0036", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("Spleen", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildUpperLargeIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildUpperLargeIntestine()", "human_phantom0037", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("UpperLargeIntestine", fTrunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLowerLargeIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLowerLargeIntestine()", "human_phantom0038", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("LowerLargeIntestine", fTrunkVolume, colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildSmallIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildSamllIntestine()", "human_phantom0039", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("SmallIntestine", fTrunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftKidney(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLeftKidney()", "human_phantom0040", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("LeftKidney", fTrunkVolume,colourName, solidVis, sensitivity);
}
void G4PhantomBuilder::BuildRightKidney(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRightKidney()", "human_phantom0041", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("RightKidney", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLeftAdrenal(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLeftAdrenal()", "human_phantom0042", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("LeftAdrenal", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildRightAdrenal(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildRightAdrenal()", "human_phantom0043", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("RightAdrenal", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildLiver(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildLiver()", "human_phantom0044", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("Liver", fTrunkVolume,colourName, solidVis, sensitivity);
}
void G4PhantomBuilder::BuildPancreas(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildPancreas()", "human_phantom0045", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("Pancreas", fTrunkVolume,colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildUrinaryBladder(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildUrinaryBladder()", "human_phantom0046", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("UrinaryBladder", fTrunkVolume, colourName, solidVis, sensitivity);
}

void G4PhantomBuilder::BuildThyroid(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (fHeadVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildThyroid()", "human_phantom0047", FatalException, "The trunk volume is missing !!!!!");

  fBody -> CreateOrgan("Thyroid", fHeadVolume, colourName,solidVis, sensitivity);
}

void G4PhantomBuilder::BuildThymus(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (fTrunkVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildThymus()", "human_phantom0048", FatalException, "The trunk volume is missing !!!!!");
  
  fBody -> CreateOrgan("Thymus",fTrunkVolume, colourName,solidVis, sensitivity);
}

G4VPhysicalVolume* G4PhantomBuilder::GetPhantom()
{
  return fMotherVolume;
}

void G4PhantomBuilder::SetMotherVolume(G4VPhysicalVolume* mother)
{
  fMotherVolume = mother;
}

void G4PhantomBuilder::SetModel(G4String modelFlag)
{
  fModel = modelFlag;

  if(fModel=="MIRD") fBody = new G4MIRDBodyFactory();
  else if(fModel=="ORNLFemale")
  {
#ifdef G4LIB_USE_GDML
    fBody = new G4ORNLFemaleBodyFactory();
#else
    G4cout << fModel << " Working with GDML only! set -DWITH_GDML_USE=ON  flag during the CMAKE configuration step" << G4endl;
#endif
  }
  else if(fModel=="ORNLMale") 
  {
#ifdef G4LIB_USE_GDML
    fBody = new G4ORNLMaleBodyFactory();
#else
    G4cout << fModel << " Working with GDML only! set -DWITH_GDML_USE=ON  flag during the CMAKE configuration step" << G4endl;
#endif
  }
}
