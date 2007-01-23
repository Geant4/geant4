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

#include "G4RunManager.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4HumanPhantomEnergyDeposit.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PVPlacement.hh"

G4HumanPhantomConstruction::G4HumanPhantomConstruction(G4HumanPhantomEnergyDeposit* energyDep):
  edepTot(energyDep)
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
  G4HumanPhantomSD* userPhantomSD = new G4HumanPhantomSD( bodypartSD, edepTot);
  SDman->AddNewDetector( userPhantomSD );
 
  G4PhantomBuilder* builder = 0;

  if(sex=="Female") builder = new G4FemaleBuilder;
  else builder = new G4MaleBuilder ;  

  builder->SetMotherVolume(ConstructWorld());

  builder->SetModel(model);

  builder->SetSex(sex);

  // the argument indicates the sensitivity of the volume
  builder->BuildHead(sensitivities["Head"], "Head", "HeadVolume", "pink", false);
  
  builder->BuildTrunk(sensitivities["Trunk"], "Trunk", "TrunkVolume", "pink", false);
   builder->BuildLeftLeg(sensitivities["LeftLeg"], "LeftLeg", "LeftLegVolume", "pink", false);
  builder->BuildRightLeg(sensitivities["RightLeg"], "RightLeg", "RightLegVolume", "pink", false);
 
  builder->BuildBrain(sensitivities["Brain"], "Brain", "BrainVolume", "pink", true); 
 
   // Remind! the elliptical cone gives problems! Intersections of volumes, 
  // wrong calculation of the volume!
  builder->BuildLeftArmBone(sensitivities["LeftArmBone"], "LeftArmBone", "LeftArmBoneVolume", "grey", true);
  builder->BuildRightArmBone(sensitivities["RightArmBone"], "RightArmBone", "RightArmBoneVolume", "grey", true);  
    

  builder->BuildLeftLegBone(sensitivities["LeftLegBone"], "LeftLegBone", "LeftLegBoneVolume", "grey", true);
  builder ->BuildRightLegBone(sensitivities["RightLegBone"], "RightLegBone", "RightLegBoneVolume", "grey", true);
  
  builder->BuildSkull(sensitivities["Skull"], "Skull", "SkullVolume", "grey", false); 
 
  builder->BuildUpperSpine(sensitivities["UpperSpine"], "UpperSpine", "UpperSpineVolume", "grey", true);
  
  builder->BuildMiddleLowerSpine(sensitivities["MiddleLowerSpine"], "MiddleLowerSpine", "MiddleLowerSpineVolume", "grey", true);
  
  builder->BuildPelvis(sensitivities["Pelvis"], "Pelvis", "PelvisVolume", "grey", true); 
  
  builder->BuildStomach(sensitivities["Stomach"],"Stomach", "StomachVolume", "yellow", true ); 
  builder->BuildUpperLargeIntestine(sensitivities["UpperLargeIntestine"],"UpperLargeIntestine", "UpperLargeIntestineVolume", "yellow", true);
  builder->BuildLowerLargeIntestine(sensitivities["LowerLargeIntestine"], "LowerLargeIntestine", "LowerLargeIntestineVolume", "yellow", true);
  if (model == "MIRD") builder->BuildRibCage(sensitivities["RibCage"], "RibCage", "RibCageVolume", "grey", true); 
  builder->BuildSpleen(sensitivities["Spleen"], "Spleen", "SpleenVolume", "brown", true);
  builder->BuildPancreas(sensitivities["Pancreas"], "Pancreas", "PancreasVolume", "lightBlue", true); 
  //builder->BuildLiver(sensitivities["Liver"]);  da fare

  builder->BuildLeftKidney(sensitivities["LeftKidney"], "LeftKidney", "LeftKidneyVolume", "green", true);
  builder->BuildRightKidney(sensitivities["RightKidney"], "RightKidney", "RightKidneyVolume", "green", true);
  builder->BuildUrinaryBladder(sensitivities["UrinaryBladder"], "UrinaryBladder", "UrinaryBladderVolume", "yellow", true);
 
  //builder->BuildHeart(sensitivities["Hearth"]);
  builder->BuildLeftLung(sensitivities["LeftLung"], "LeftLung", "LeftLungVolume", "blue", true);
  builder->BuildRightLung(sensitivities["RightLung"], "RightLung", "RightLungVolume", "blue", true);
  //builder->BuildThyroid(sensitivities["Thyroid"]); 
  
  if(sex=="Female"){
    builder->BuildLeftOvary(sensitivities["LeftOvary"], "LeftOvary", "LeftOvaryVolume", "purple", true);
    builder->BuildRightOvary(sensitivities["RightOvary"], "RightOvary", "RightOvaryVolume", "purple", true);
    builder->BuildUterus(sensitivities["Uterus"], "Uterus", "UterusVolume", "purple", true);
    builder->BuildLeftBreast(sensitivities["LeftBreast"], "LeftBreast", "LeftBreastVolume", "purple", true); 
    builder->BuildRightBreast(sensitivities["RightBreast"], "RightBreast", "RightBreastVolume", "purple", true);
   //if(model == "MIRD") builder->BuildParameterisedBreast(true);
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

void G4HumanPhantomConstruction::CleanPhantom()
{
  //delete mother;
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
}
