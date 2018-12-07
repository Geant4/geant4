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
// Previous authors: G. Guerrieri, S. Guatelli, and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//

#include <map>

#include "globals.hh"

#include "G4HumanPhantomConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4HumanPhantomSD.hh"
#include "G4SDManager.hh"

//#include "G4VBodyFactory.hh"
//#include "G4MIRDBodyFactory.hh"
//#include "G4ORNLBodyFactory.hh"

#include "G4PhantomBuilder.hh"
#include "G4FemaleBuilder.hh"
#include "G4MaleBuilder.hh"
#include "G4PhantomHeadBuilder.hh"
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
  
  

  G4BasePhantomBuilder*  builder = 0;

  if (model == "MIRDHead" || model == "ORNLHead")
    { 
      G4cout << "HeadBuilder instantiated" << G4endl;
      builder = new G4PhantomHeadBuilder;
      if (model ==  "MIRDHead") builder->SetModel("MIRD");
      else if (model ==  "ORNLHead") builder->SetModel("ORNLMale");
    }
  else
    {  
      if (sex =="Female") 
	{ 
	  if (model == "MIX") builder = new G4CustomFemaleBuilder;
	  else {builder = new G4FemaleBuilder;}
	  builder->SetModel(model);
	  G4cout <<model << " "<< sex << G4endl;
	}
      else if (sex == "Male") 
	{
	  builder = new G4MaleBuilder;
	  builder->SetModel(model);
          if (model == "MIX") 
	    { 
	      G4cout<< "Custom Male is not available!!! MIRD model is selected !" 
		    << G4endl;
	      model = "MIRD";  
	      builder->SetModel(model);}
	}
    }
  
  builder->SetMotherVolume(ConstructWorld());
  
  // the argument indicates the sensitivity of the volume
  
  builder->BuildHead("black", false, sensitivities["Head"]);
  builder->BuildSkull("orange", false,sensitivities["Skull"]); 
  builder->BuildBrain("yellow", true,sensitivities["Brain"]); 

  if (model != "MIRDHead" && model != "ORNLHead")
    { 
      //  builder->SetModel(model);
      builder->BuildTrunk("yellow", false, sensitivities["Trunk"]);
      
      builder->BuildLeftLeg("yellow", false,sensitivities["LeftLeg"]);
      builder->BuildRightLeg("yellow", false,sensitivities["RightLeg"]);
      
      builder->BuildLeftArmBone("grey", true,sensitivities["LeftArmBone"]);
      builder->BuildRightArmBone("grey", true, sensitivities["RightArmBone"]);  
    
      builder->BuildLeftLegBone("grey", true,sensitivities["LeftLegBone"]);
      builder ->BuildRightLegBone("grey", true,sensitivities["RightLegBone"]);
 
      builder->BuildUpperSpine("yellow", true,sensitivities["UpperSpine"]); 
      
      if (model == "MIRD" || model == "MIX") 
	{
	  builder->BuildLeftScapula("grey", true, sensitivities["LeftScapula"]); 
	  builder->BuildRightScapula("grey", true, sensitivities["RightScapula"]);
	  builder->BuildLeftAdrenal("yellow", true, sensitivities["LeftAdrenal"]);
	  builder->BuildRightAdrenal("yellow", true, sensitivities["RightAdrenal"]);
	  builder->BuildThymus("orange", true,sensitivities["Thymus"]); 
	  builder->BuildLeftClavicle("grey", true,sensitivities["LeftClavicle"]);
	  builder->BuildRightClavicle("grey", true,sensitivities["RightClavicle"]);
	  builder->BuildSmallIntestine("orange", true,sensitivities["SmallIntestine"]);
	  builder->BuildRibCage("grey", true,sensitivities["RibCage"]); 
	}
  
      builder->BuildMiddleLowerSpine("yellow", true,sensitivities["MiddleLowerSpine"]);
  
      builder->BuildPelvis("grey", true,sensitivities["Pelvis"]); 
  
      builder->BuildStomach("orange", true,sensitivities["Stomach"]); 
      builder->BuildUpperLargeIntestine("lightBlue", true,sensitivities["UpperLargeIntestine"]);
      builder->BuildLowerLargeIntestine("lightBlue", true,sensitivities["LowerLargeIntestine"]);
         
      builder->BuildSpleen("green", true,sensitivities["Spleen"]);
      builder->BuildPancreas("purple", true,sensitivities["Pancreas"]); 
      //builder->BuildLiver("orange", true,sensitivities["Liver"]); 

      builder->BuildLeftKidney("green", true,sensitivities["LeftKidney"]);
      builder->BuildRightKidney("green", true,sensitivities["RightKidney"]);
      builder->BuildUrinaryBladder("green", true,sensitivities["UrinaryBladder"]);
 
      //builder->BuildHeart("red", true,sensitivities["Hearth"]);// to do MIRD
     // builder->BuildLeftLung("blue", true,sensitivities["LeftLung"]);
      //builder->BuildRightLung("blue", true,sensitivities["RightLung"]);
     // builder->BuildThyroid("orange", true,sensitivities["Thyroid"]); 

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
	    builder->BuildVoxelRightBreast("purple", false, sensitivities["RightBreast"]);  
	  } 
      }
      
      if(sex=="Male"){
	
	if (model == "MIRD"){ 
	  builder -> BuildMaleGenitalia("yellow",false,sensitivities["MaleGenitalia"]);
	  builder -> BuildLeftTeste("purple",true,sensitivities["LeftTeste"]);
	  builder -> BuildRightTeste("purple",true,sensitivities["RightTeste"]);
	}
	else G4cout <<  "ORNL does not have model for male genitalia and testes yet" << G4endl;
      }
      
    }
  G4VPhysicalVolume* result=builder->GetPhantom(); 
  delete builder;
  return result; 
}

void  G4HumanPhantomConstruction::SetBodyPartSensitivity(G4String, G4bool)
{
  G4cout << "This method is not currently working !!!!" << G4endl;
}

G4VPhysicalVolume* G4HumanPhantomConstruction::ConstructWorld()
{
  G4Material* air = material -> GetMaterial("Air");

  // World Volume
  // G4double worldSize = 1.*m ;
  G4double worldSize = 1.5 *m ;
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
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
 
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
  if (model == "MIRDHead")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }

  if (model == "ORNLHead")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }
}

void G4HumanPhantomConstruction::ConstructSDandField()
{
   G4HumanPhantomSD* SD = new G4HumanPhantomSD("SD", "HumanPhantomCollection");
   G4SDManager::GetSDMpointer()->AddNewDetector(SD);
if (model != "ORNLMale" && model != "ORNLFemale" && model!= "ORNLHead")  
{
  SetSensitiveDetector("logicalHead",SD);
  SetSensitiveDetector("logicalSkull",SD);
  SetSensitiveDetector("logicalBrain",SD);
  if (model != "MIRDHead")
    { 
      SetSensitiveDetector("logicalTrunk",SD);
      SetSensitiveDetector("logicalLeftLeg",SD); 
      SetSensitiveDetector("logicalRightLeg",SD);
      SetSensitiveDetector("logicalLeftArmBone",SD); 
      SetSensitiveDetector("logicalRightArmBone",SD);
      SetSensitiveDetector("logicalLeftLegBone",SD); 
      SetSensitiveDetector("logicalRightLegBone",SD);
      SetSensitiveDetector("logicalUpperSpine",SD);
      SetSensitiveDetector("logicalLeftScapula",SD);
      SetSensitiveDetector("logicalRightScapula",SD);
      SetSensitiveDetector("logicalLeftAdrenal",SD);
      SetSensitiveDetector("logicalRightAdrenal",SD);      SetSensitiveDetector("logicalThymus",SD);      SetSensitiveDetector("logicalLeftClavicle",SD);
      SetSensitiveDetector("logicalRightClavicle",SD);
      SetSensitiveDetector("logicalSmallIntestine",SD); 
      SetSensitiveDetector("logicalRibCage",SD);       SetSensitiveDetector("logicalMiddleLowerSpine",SD); 
      SetSensitiveDetector("logicalStomach",SD); 
      SetSensitiveDetector("logicalUpperLargeIntestine",SD);
      SetSensitiveDetector("logicalLowerLargeIntestine",SD);
      SetSensitiveDetector("logicalSpleen",SD);
      SetSensitiveDetector("logicalPancreas",SD);
      SetSensitiveDetector("logicalLeftKidney",SD);
      SetSensitiveDetector("logicalRightKidney",SD);       
      SetSensitiveDetector("logicalUrinaryBladder",SD);

      if(sex=="Female"){

	SetSensitiveDetector("logicalLeftOvary",SD);
        SetSensitiveDetector("logicalRightOvary",SD); 
        SetSensitiveDetector("logicalUterus",SD);
        SetSensitiveDetector("logicalLeftBreast",SD);
        SetSensitiveDetector("logicalRightBreast",SD); 
	}

      if(sex=="Male"){
	
	
	  SetSensitiveDetector("logicalMaleGenitalia",SD);
          SetSensitiveDetector("logicalLeftTeste",SD);
	  SetSensitiveDetector("logicalRightTeste",SD);
	}
	
      }
  }else 
 { 
  SetSensitiveDetector("HeadVolume",SD);
  SetSensitiveDetector("SkullVolume",SD);
  SetSensitiveDetector("BrainVolume",SD);
  G4cout << "Work in progress!!!! " << G4endl;
  G4cout <<"ORNL model!!!! Head is sensitive only!!!" << G4endl;
} 

   
}
