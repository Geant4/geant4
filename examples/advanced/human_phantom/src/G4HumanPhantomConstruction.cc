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
#include "globals.hh"
#include <map>

#include "G4HumanPhantomConstruction.hh"

#include "G4HumanPhantomSD.hh"
#include "G4SDManager.hh"

#include "G4VBodyFactory.hh"
#include "G4MIRDBodyFactory.hh"
#include "G4ORNLBodyFactory.hh"

#include "G4FemaleBuilder.hh"
#include "G4MaleBuilder.hh"

#include "G4RunManager.hh"


G4HumanPhantomConstruction::G4HumanPhantomConstruction()
{
 messenger = new G4HumanPhantomMessenger(this);
}

G4HumanPhantomConstruction::~G4HumanPhantomConstruction()
{
 delete messenger;
 delete userPhantomSD;
}


G4VPhysicalVolume* G4HumanPhantomConstruction::Construct()
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String bodypartSD = "BodyPartSD";
  G4HumanPhantomSD* userPhantomSD = new G4HumanPhantomSD( bodypartSD );
  SDman->AddNewDetector( userPhantomSD );

  G4FemaleBuilder builder; 
  //G4MaleBuilder builder; 

  builder.BuildWorld();

  builder.SetModel(model);

  builder.SetSex(sex);

  builder.BuildHead(sensitivities["Head"]);
  builder.BuildTrunk(sensitivities["Trunk"]);
  builder.BuildLegs(sensitivities["Legs"]);

  builder.BuildBreast(sensitivities["Breast"]);

  builder.BuildBrain(sensitivities["Brain"]);

  builder.BuildArmBone(sensitivities["ArmBone"]);
  builder.BuildLegBone(sensitivities["LegBone"]);
  builder.BuildSkull(sensitivities["Skull"]);
  builder.BuildUpperSpine(sensitivities["UpperSpine"]);
  builder.BuildMiddleLowerSpine(sensitivities["MiddleLowerSpine"]);
  builder.BuildPelvis(sensitivities["Pelvis"]);

  builder.BuildStomach(sensitivities["Stomach"]);
  builder.BuildUpperLargeIntestine(sensitivities["UpperLargeIntestine"]);
  builder.BuildLowerLargeIntestine(sensitivities["LowerLargeIntestine"]);

  builder.BuildSpleen(sensitivities["Spleen"]);
  builder.BuildPancreas(sensitivities["Pancreas"]);
  builder.BuildLiver(sensitivities["Liver"]);
  builder.BuildKidney(sensitivities["Kidney"]);
  builder.BuildUrinaryBladder(sensitivities["UrinaryBladder"]);

  builder.BuildHeart(sensitivities["Heart"]);

  builder.BuildLung(sensitivities["Lung"]);

  builder.BuildThyroid(sensitivities["Thyroid"]);

  builder.BuildOvary(sensitivities["Ovary"]);
  builder.BuildUterus(sensitivities["Uterus"]);

  return builder.GetPhantom(); 
}

void  G4HumanPhantomConstruction::SetBodyPartSensitivity(G4String bodyPartName, G4bool bodyPartSensitivity)
{
  sensitivities[bodyPartName] = bodyPartSensitivity;
  if(bodyPartSensitivity==true) 
  G4cout << " >>> " << bodyPartName << " added as sensitive volume." << G4endl;
}

void G4HumanPhantomConstruction::CleanPhantom()
{
 delete mother;
}

void G4HumanPhantomConstruction::UpdatePhantom()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

void G4HumanPhantomConstruction::SetPhantomSex(G4String newSex)
{
  sex=newSex;

  if (sex == "Male")
    {
      G4cout << ">> Male Phantom will be built." << G4endl;
    }
  if (sex == "Female")
    {
      G4cout << ">> Female Phantom will be built." << G4endl;
    }
  if ((sex != "Female") && (sex != "Male"))
    G4cout << sex << " can not be defined!" << G4endl;
}

void G4HumanPhantomConstruction::SetPhantomModel(G4String newModel)
{
  model = newModel;

  if (model == "MIRD")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }
  if (model == "ORNL")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }
  if (model == "MIX")
    {
      G4cout<<" >> Phantom with mixed body parts will be built."<<G4endl;
    }
}
