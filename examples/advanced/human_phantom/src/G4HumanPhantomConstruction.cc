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

#include "globals.hh"
#include <map>

#include "G4HumanPhantomConstruction.hh"

#include "G4HumanPhantomSD.hh"
#include "G4SDManager.hh"

//#include "G4VBodyFactory.hh"
//#include "G4MIRDBodyFactory.hh"
//#include "G4ORNLBodyFactory.hh"

#include "G4PhantomBuilder.hh"
#include "G4FemaleBuilder.hh"
#include "G4MaleBuilder.hh"
#include "G4CustomFemaleBuilder.hh"
#include "G4RunManager.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PVPlacement.hh"

G4HumanPhantomConstruction::G4HumanPhantomConstruction()
{
 messenger = new G4HumanPhantomMessenger(this);
 material = new G4HumanPhantomMaterial();
}

G4HumanPhantomConstruction::~G4HumanPhantomConstruction()
{
 delete material;
 delete messenger;
}

G4VPhysicalVolume* G4HumanPhantomConstruction::Construct()
{
   material -> DefineMaterials();
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String bodypartSD = "BodyPartSD";
  G4HumanPhantomSD* userPhantomSD = new G4HumanPhantomSD(bodypartSD);
  SDman->AddNewDetector(userPhantomSD);
 
  //G4PhantomBuilder::SetModel(model);
  G4PhantomBuilder* builder = 0;
  
  if (sex =="Female") 
    { 
      if (model == "MIX") builder = new G4CustomFemaleBuilder;
         else builder = new G4FemaleBuilder;
    }
  else if (sex == "Male") builder = new G4MaleBuilder ;
  else if (sex == "Male" && model == "MIX") 
    { 
      G4cout<< "Custom Male is not available!!! MIRD model is selected !" << G4endl;
      model = "MIRD";  
      builder = new G4MaleBuilder ;} 

  builder->SetMotherVolume(ConstructWorld());

  builder->SetModel(model);

  // the argument indicates the sensitivity of the volume

  builder->BuildHead("pink", false, sensitivities["Head"]);
  
  builder->BuildTrunk("pink", false, sensitivities["Trunk"]);
  
  builder->BuildLeftLeg("pink", false,sensitivities["LeftLeg"]);
  builder->BuildRightLeg("pink", false,sensitivities["RightLeg"]);
  builder->BuildBrain("pink", true,sensitivities["Brain"]); 
 
  builder->BuildLeftArmBone("grey", true,sensitivities["LeftArmBone"]);
  builder->BuildRightArmBone("grey", true, sensitivities["RightArmBone"]);  
    
  builder->BuildLeftLegBone("grey", true,sensitivities["LeftLegBone"]);
  builder ->BuildRightLegBone("grey", true,sensitivities["RightLegBone"]);
  
  builder->BuildSkull("grey", false,sensitivities["Skull"]); 
 
  builder->BuildUpperSpine("grey", true,sensitivities["UpperSpine"]);
  
  builder->BuildMiddleLowerSpine("grey", true,sensitivities["MiddleLowerSpine"]);
  
  builder->BuildPelvis("grey", true,sensitivities["Pelvis"]); 
  
  builder->BuildStomach("yellow", true,sensitivities["Stomach"]); 
  builder->BuildUpperLargeIntestine("yellow", true,sensitivities["UpperLargeIntestine"]);
  builder->BuildLowerLargeIntestine("yellow", true,sensitivities["LowerLargeIntestine"]);
  if (model == "MIRD") builder->BuildRibCage("grey", true,sensitivities["RibCage"]); 
  builder->BuildSpleen("brown", true,sensitivities["Spleen"]);
  builder->BuildPancreas("lightBlue", true,sensitivities["Pancreas"]); 
  //builder->BuildLiver(sensitivities["Liver"]);  da fare

  builder->BuildLeftKidney("green", true,sensitivities["LeftKidney"]);
  builder->BuildRightKidney("green", true,sensitivities["RightKidney"]);
  builder->BuildUrinaryBladder("yellow", true,sensitivities["UrinaryBladder"]);
 
  //builder->BuildHeart(sensitivities["Hearth"]);
  builder->BuildLeftLung("blue", true,sensitivities["LeftLung"]);
  builder->BuildRightLung("blue", true,sensitivities["RightLung"]);
  //builder->BuildThyroid(sensitivities["Thyroid"]); 
  
  if(sex=="Female"){

    builder->BuildLeftOvary("purple", true,sensitivities["LeftOvary"]);
    builder->BuildRightOvary("purple", true,sensitivities["RightOvary"]);
    builder->BuildUterus("purple", true,sensitivities["Uterus"]);

    if (model == "ORNLFemale" || model == "MIRD")
      {
       builder->BuildLeftBreast("purple", true,sensitivities["LeftBreast"]); 
       builder->BuildRightBreast("purple", true,sensitivities["RightBreast"]);
      }
    else if (model == "MIX")
      {
       builder->BuildVoxelLeftBreast("purple",false, sensitivities["LeftBreast"]); 
       //builder->BuildParameterisedRightBreast("purple", false, sensitivities["RightBreast"]);  
      } 
 }
  /*  
  if(sex=="Male"){
    //    builder->BuildMaleGenitalia(sensitivities["MaleGenitalia"]);
    // builder->BuildTestes(sensitivities["Testes"]);
  }
  */
  return builder->GetPhantom(); 
  delete builder;
}

void  G4HumanPhantomConstruction::SetBodyPartSensitivity(G4String bodyPartName, G4bool bodyPartSensitivity)
{
   sensitivities[bodyPartName] = bodyPartSensitivity;
  if(bodyPartSensitivity==true) 
  G4cout << " >>> " << bodyPartName << " added as sensitive volume." << G4endl;
}

G4VPhysicalVolume* G4HumanPhantomConstruction::ConstructWorld()
{
  G4Material* air = material -> GetMaterial("Air");

// World Volume
  G4double worldSize = 1.*m ;

  G4Box* world = new G4Box("world", worldSize, worldSize, worldSize);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(world, 
						    air, 
						    "logicalWorld", 0, 0,0);

  G4VPhysicalVolume* motherVolume = new G4PVPlacement(0,G4ThreeVector(),
						"physicalWorld",
						logicWorld,
						0,
						false,
						0);

  // Visualization Attributes
 G4VisAttributes* WorldVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
  
  WorldVisAtt->SetForceSolid(false);
  logicWorld->SetVisAttributes(WorldVisAtt);
 
 return motherVolume;
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
  if (model == "ORNLFemale")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }

if (model == "ORNLMale")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }

if (model == "MIX")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }
}
