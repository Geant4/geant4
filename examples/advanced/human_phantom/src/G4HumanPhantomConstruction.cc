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
// This is the Director of the Builder: 
//  define sex at line 62!
//  default is Female.
//
//*********************************************************************
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

G4HumanPhantomConstruction::G4HumanPhantomConstruction()
{
 messenger = new G4HumanPhantomMessenger(this);
 material = new G4HumanPhantomMaterial();
}

G4HumanPhantomConstruction::~G4HumanPhantomConstruction()
{
 delete material;
 delete messenger;
 delete userPhantomSD;
}


G4VPhysicalVolume* G4HumanPhantomConstruction::Construct()
{
  material -> DefineMaterials();

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String bodypartSD = "BodyPartSD";
  G4HumanPhantomSD* userPhantomSD = new G4HumanPhantomSD( bodypartSD );
  SDman->AddNewDetector( userPhantomSD );

  G4String model = "MIRD";
   // Define sex of phantom
  G4String sex = "Female";

  G4FemaleBuilder builder ;  

  builder.BuildWorld();

  builder.SetModel(model);

  builder.SetSex(sex);

  // the argument indicates the sensitivity of the volume
  builder.BuildHead(true);
  builder.BuildTrunk(true);
  builder.BuildLegs(true);

  
  builder.BuildBrain(true);
  builder.BuildArmBone(true);
  builder.BuildLegBone(true);

  
  builder.BuildSkull(true);
  
  builder.BuildUpperSpine(true);

  builder.BuildMiddleLowerSpine(true);

  builder.BuildPelvis(true);
  builder.BuildStomach(true);
 
  
  builder.BuildUpperLargeIntestine(true);
  builder.BuildLowerLargeIntestine(true);
  
  builder.BuildSpleen(true);
  builder.BuildPancreas(true); 

  builder.BuildLiver(true); // da controllare
  
  builder.BuildKidney(true); 
 
  builder.BuildUrinaryBladder(true);

  builder.BuildHeart(true); // controllare 
 
  builder.BuildLung(true);

  builder.BuildThyroid(true); 

  
  if(sex=="Female"){
    builder.BuildOvary(true);
    builder.BuildUterus(true);
    builder.BuildBreast(true);
    builder.BuildParameterisedBreast(true);
  }
  /*
  if(sex=="Male"){
    // builder.BuildMaleGenitalia(sensitivities["MaleGenitalia"]);
    // builder.BuildTestes(sensitivities["Testes"]);
  }
  */
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
