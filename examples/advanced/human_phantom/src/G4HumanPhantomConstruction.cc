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

  builder->BuildWorld();

  builder->SetModel(model);

  builder->SetSex(sex);

  // the argument indicates the sensitivity of the volume
  builder->BuildHead(sensitivities["Head"]);

  builder->BuildTrunk(sensitivities["Trunk"]);
  builder->BuildLegs(sensitivities["Legs"]);
  builder->BuildBrain(sensitivities["Brain"]);

  builder->BuildArmBone(sensitivities["ArmBone"]);
  builder->BuildLegBone(sensitivities["LegBone"]);
  builder->BuildSkull(sensitivities["Skull"]);
  builder->BuildUpperSpine(sensitivities["UpperSpine"]);
  builder->BuildMiddleLowerSpine(sensitivities["MiddleLowerSpine"]);

  builder->BuildPelvis(sensitivities["Pelvis"]);
  builder->BuildStomach(sensitivities["Stomach"]);
  builder->BuildUpperLargeIntestine(sensitivities["UpperLargeIntestine"]);
  builder->BuildLowerLargeIntestine(sensitivities["LowerLargeIntestine"]);

  builder->BuildSpleen(sensitivities["Spleen"]);
  builder->BuildPancreas(sensitivities["Pancreas"]); 
  builder->BuildLiver(sensitivities["Liver"]);  

  builder->BuildKidney(sensitivities["Kidney"]); 
  builder->BuildUrinaryBladder(sensitivities["UrinaryBladder"]);

  builder->BuildHeart(sensitivities["Hearth"]);
  builder->BuildLung(sensitivities["Lung"]);
  builder->BuildThyroid(sensitivities["Thyroid"]); 

  if(sex=="Female"){
    builder->BuildOvary(true);
    builder->BuildUterus(true);
    builder->BuildBreast(true);
    if(model == "MIRD") builder->BuildParameterisedBreast(true);
  }
 
  if(sex=="Male"){
     builder->BuildMaleGenitalia(sensitivities["MaleGenitalia"]);
     builder->BuildTestes(sensitivities["Testes"]);
  }
  
  return builder->GetPhantom(); 
  delete builder;
}

void  G4HumanPhantomConstruction::SetBodyPartSensitivity(G4String bodyPartName, G4bool bodyPartSensitivity)
{
   sensitivities[bodyPartName] = bodyPartSensitivity;
  if(bodyPartSensitivity==true) 
  G4cout << " >>> " << bodyPartName << " added as sensitive volume." << G4endl;
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
