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
  fMessenger = new G4HumanPhantomMessenger(this);
  fMaterial = new G4HumanPhantomMaterial();
}

G4HumanPhantomConstruction::~G4HumanPhantomConstruction()
{
  delete fMaterial;
  delete fMessenger;
}

G4VPhysicalVolume* G4HumanPhantomConstruction::Construct()
{
  fMaterial -> DefineMaterials();
  
  G4BasePhantomBuilder*  builder = nullptr;

  if (fModel == "MIRDHead" || fModel == "ORNLHead")
    { 
      G4cout << "HeadBuilder instantiated" << G4endl;
      builder = new G4PhantomHeadBuilder;
      if (fModel ==  "MIRDHead") builder->SetModel("MIRD");
      else if (fModel ==  "ORNLHead") builder->SetModel("ORNLMale");
    }
  else
    {  
      if (fSex =="Female") 
	{ 
	  builder = new G4FemaleBuilder;
	  builder->SetModel(fModel);
	  G4cout << fModel << " "<< fSex << G4endl;
	}
      else if (fSex == "Male") 
	{
	  builder = new G4MaleBuilder;
	  builder->SetModel(fModel);
	}
    }
  
  builder->SetMotherVolume(ConstructWorld());
  
  // the argument indicates the sensitivity of the volume
  
  builder->BuildHead("black", false, fSensitivities["Head"]);
  builder->BuildSkull("orange", false, fSensitivities["Skull"]); 
  builder->BuildBrain("yellow", true, fSensitivities["Brain"]); 

  if (fModel != "MIRDHead" && fModel != "ORNLHead")
    { 
      //  builder->SetModel(model);
      builder->BuildTrunk("yellow", false, fSensitivities["Trunk"]);
      
      builder->BuildLeftLeg("yellow", false, fSensitivities["LeftLeg"]);
      builder->BuildRightLeg("yellow", false, fSensitivities["RightLeg"]);
      
      builder->BuildLeftArmBone("grey", true, fSensitivities["LeftArmBone"]);
      builder->BuildRightArmBone("grey", true, fSensitivities["RightArmBone"]);  
    
      builder->BuildLeftLegBone("grey", true, fSensitivities["LeftLegBone"]);
      builder ->BuildRightLegBone("grey", true, fSensitivities["RightLegBone"]);
 
      builder->BuildUpperSpine("yellow", true, fSensitivities["UpperSpine"]); 
      
      if (fModel == "MIRD") 
	{
	  builder->BuildLeftScapula("grey", true, fSensitivities["LeftScapula"]); 
	  builder->BuildRightScapula("grey", true, fSensitivities["RightScapula"]);
	  builder->BuildLeftAdrenal("yellow", true, fSensitivities["LeftAdrenal"]);
	  builder->BuildRightAdrenal("yellow", true, fSensitivities["RightAdrenal"]);
	  builder->BuildThymus("orange", true, fSensitivities["Thymus"]); 
	  builder->BuildLeftClavicle("grey", true, fSensitivities["LeftClavicle"]);
	  builder->BuildRightClavicle("grey", true, fSensitivities["RightClavicle"]);
	  builder->BuildSmallIntestine("orange", true, fSensitivities["SmallIntestine"]);
	  builder->BuildRibCage("grey", true, fSensitivities["RibCage"]); 
	}
  
      builder->BuildMiddleLowerSpine("yellow", true, fSensitivities["MiddleLowerSpine"]);
  
      builder->BuildPelvis("grey", true, fSensitivities["Pelvis"]); 
  
      builder->BuildStomach("orange", true, fSensitivities["Stomach"]); 
      builder->BuildUpperLargeIntestine("lightBlue", true, fSensitivities["UpperLargeIntestine"]);
      builder->BuildLowerLargeIntestine("lightBlue", true, fSensitivities["LowerLargeIntestine"]);
      builder->BuildSpleen("green", true, fSensitivities["Spleen"]);
      builder->BuildPancreas("purple", true, fSensitivities["Pancreas"]); 
      //builder->BuildLiver("orange", true,sensitivities["Liver"]); 

      builder->BuildLeftKidney("green", true, fSensitivities["LeftKidney"]);
      builder->BuildRightKidney("green", true, fSensitivities["RightKidney"]);
      builder->BuildUrinaryBladder("green", true, fSensitivities["UrinaryBladder"]);
 
      //builder->BuildHeart("red", true,fSensitivities["Hearth"]);// to do MIRD
     // builder->BuildLeftLung("blue", true, fSensitivities["LeftLung"]);
      //builder->BuildRightLung("blue", true, fSensitivities["RightLung"]);
     // builder->BuildThyroid("orange", true, fSensitivities["Thyroid"]); 

      if(fSex=="Female"){

	builder->BuildLeftOvary("purple", true, fSensitivities["LeftOvary"]);
	builder->BuildRightOvary("purple", true, fSensitivities["RightOvary"]);
	builder->BuildUterus("purple", true, fSensitivities["Uterus"]);

	if (fModel == "ORNLFemale" || fModel == "MIRD")
	  {
	    builder->BuildLeftBreast("purple", true,fSensitivities["LeftBreast"]); 
	    builder->BuildRightBreast("purple", true,fSensitivities["RightBreast"]);
	  }
      }
      
      if(fSex=="Male"){
	
	if (fModel == "MIRD"){ 
	  builder -> BuildMaleGenitalia("yellow",false, fSensitivities["MaleGenitalia"]);
	  builder -> BuildLeftTeste("purple",true, fSensitivities["LeftTeste"]);
	  builder -> BuildRightTeste("purple",true, fSensitivities["RightTeste"]);
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
  G4Material* air = fMaterial -> GetMaterial("Air");

  // World Volume
  // G4double worldSize = 1.*m ;
  G4double worldSize = 1.5 *m ;
  G4Box* world = new G4Box("world", worldSize, worldSize, worldSize);

  auto* logicWorld = new G4LogicalVolume(world, 
				          air, 
				          "logicalWorld", nullptr, nullptr, nullptr);

  G4VPhysicalVolume* motherVolume = new G4PVPlacement(nullptr,G4ThreeVector(),
						      "physicalWorld",
						      logicWorld,
						      nullptr,
						      false,
						      0);

  // Visualization Attributes
  auto* WorldVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
    
  WorldVisAtt->SetForceSolid(false);
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
 
  return motherVolume;
}

void G4HumanPhantomConstruction::SetPhantomSex(G4String newSex)
{
  fSex=newSex;

  if (fSex == "Male")
    {
      G4cout << ">> Male Phantom will be built." << G4endl;
    }
  if (fSex == "Female")
    {
      G4cout << ">> Female Phantom will be built." << G4endl;
    }
  if ((fSex != "Female") && (fSex != "Male"))
    G4cout << fSex << " can not be defined!" << G4endl;
}

void G4HumanPhantomConstruction::SetPhantomModel(G4String newModel)
{
  fModel = newModel;

  if (fModel == "MIRD")
    {
      G4cout<<" >> Phantom " << fModel << " will be built."<<G4endl;
    }
  if (fModel == "ORNLFemale")
    {
      G4cout<<" >> Phantom " << fModel << " will be built."<<G4endl;
    }

  if (fModel == "ORNLMale")
    {
      G4cout<<" >> Phantom " << fModel << " will be built."<<G4endl;
    }

  if (fModel == "MIRDHead")
    {
      G4cout<<" >> Phantom " << fModel << " will be built."<<G4endl;
    }

  if (fModel == "ORNLHead")
    {
      G4cout<<" >> Phantom " << fModel << " will be built."<<G4endl;
    }
}

void G4HumanPhantomConstruction::ConstructSDandField()
{
   auto* SD = new G4HumanPhantomSD("SD", "HumanPhantomCollection");
   G4SDManager::GetSDMpointer()->AddNewDetector(SD);

if (fModel != "ORNLMale" && fModel != "ORNLFemale" && fModel!= "ORNLHead")  
{
  SetSensitiveDetector("logicalHead",SD);
  SetSensitiveDetector("logicalSkull",SD);
  SetSensitiveDetector("logicalBrain",SD);
  if (fModel != "MIRDHead")
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
      SetSensitiveDetector("logicalRightAdrenal",SD);      
      SetSensitiveDetector("logicalThymus",SD);      
      SetSensitiveDetector("logicalLeftClavicle",SD);
      SetSensitiveDetector("logicalRightClavicle",SD);
      SetSensitiveDetector("logicalSmallIntestine",SD); 
      SetSensitiveDetector("logicalRibCage",SD);       
      SetSensitiveDetector("logicalMiddleLowerSpine",SD); 
      SetSensitiveDetector("logicalStomach",SD); 
      SetSensitiveDetector("logicalUpperLargeIntestine",SD);
      SetSensitiveDetector("logicalLowerLargeIntestine",SD);
      SetSensitiveDetector("logicalSpleen",SD);
      SetSensitiveDetector("logicalPancreas",SD);
      SetSensitiveDetector("logicalLeftKidney",SD);
      SetSensitiveDetector("logicalRightKidney",SD);       
      SetSensitiveDetector("logicalUrinaryBladder",SD);

      if(fSex=="Female"){
	SetSensitiveDetector("logicalLeftOvary",SD);
        SetSensitiveDetector("logicalRightOvary",SD); 
        SetSensitiveDetector("logicalUterus",SD);
        SetSensitiveDetector("logicalLeftBreast",SD);
        SetSensitiveDetector("logicalRightBreast",SD); 
	}
      else if(fSex=="Male"){
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
  G4cout <<"ORNL model!!!! Head is sensitive only!!!" << G4endl;
}   
}
