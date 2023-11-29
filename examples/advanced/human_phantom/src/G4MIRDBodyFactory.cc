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
//
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//
//
#include "G4MIRDBodyFactory.hh"
#include "G4MIRDStomach.hh"
#include "G4MIRDSmallIntestine.hh"
#include "G4MIRDUpperLargeIntestine.hh"
#include "G4MIRDLowerLargeIntestine.hh"
#include "G4MIRDLeftKidney.hh"
#include "G4MIRDRightKidney.hh"
#include "G4MIRDLiver.hh"
#include "G4MIRDPancreas.hh"
#include "G4MIRDSpleen.hh"
#include "G4MIRDUrinaryBladder.hh"
#include "G4MIRDLeftLung.hh"
#include "G4MIRDRightLung.hh"
#include "G4MIRDHeart.hh"
#include "G4MIRDBrain.hh"
#include "G4MIRDHead.hh"
#include "G4MIRDTrunk.hh"
#include "G4MIRDMaleGenitalia.hh"
#include "G4MIRDLeftLeg.hh"
#include "G4MIRDRightLeg.hh"
#include "G4MIRDThyroid.hh"
#include "G4MIRDThymus.hh"
#include "G4MIRDUterus.hh"
#include "G4MIRDLeftBreast.hh"
#include "G4MIRDRightBreast.hh"
#include "G4MIRDRightOvary.hh"
#include "G4MIRDLeftOvary.hh"
#include "G4MIRDUpperSpine.hh"
#include "G4MIRDMiddleLowerSpine.hh"
#include "G4MIRDLeftLegBone.hh"
#include "G4MIRDRightLegBone.hh"
#include "G4MIRDLeftClavicle.hh"
#include "G4MIRDRightClavicle.hh"
#include "G4MIRDLeftArmBone.hh"
#include "G4MIRDRightArmBone.hh"
#include "G4MIRDSkull.hh"
#include "G4MIRDRibCage.hh"
#include "G4MIRDPelvis.hh"
#include "G4MIRDLeftTeste.hh"
#include "G4MIRDRightTeste.hh"
#include "G4MIRDLeftScapula.hh"
#include "G4MIRDRightScapula.hh"
#include "G4MIRDLeftAdrenal.hh"
#include "G4MIRDRightAdrenal.hh"

G4MIRDBodyFactory::G4MIRDBodyFactory()
{
  // Map with name of the fOrgan and pointer to the MIRDOrgan class
  fOrgan["Head"] = new G4MIRDHead();
  fOrgan["Trunk"] = new G4MIRDTrunk(); 
  fOrgan["LeftLeg"] = new G4MIRDLeftLeg();
  fOrgan["RightLeg"] = new G4MIRDRightLeg();
  fOrgan["Skull"] = new G4MIRDSkull();
  fOrgan["LeftArmBone"] = new G4MIRDLeftArmBone();
  fOrgan["RightArmBone"] = new G4MIRDRightArmBone();
  fOrgan["UpperSpine"] = new G4MIRDUpperSpine();
  fOrgan["MiddleLowerSpine"] = new G4MIRDMiddleLowerSpine();
  fOrgan["Pelvis"]= new G4MIRDPelvis();
  fOrgan["RibCage"] = new G4MIRDRibCage(); 
  fOrgan["LeftClavicle"]= new G4MIRDLeftClavicle();
  fOrgan["RightClavicle"] = new G4MIRDRightClavicle();
  fOrgan["LeftLegBone"] = new G4MIRDLeftLegBone();
  fOrgan["RightLegBone"] = new G4MIRDRightLegBone();
  fOrgan["LeftScapula"]= new G4MIRDLeftScapula(); 
  fOrgan["RightScapula"]= new G4MIRDRightScapula(); 
  fOrgan["Heart"] = new G4MIRDHeart();
  fOrgan["Thyroid"] = new G4MIRDThyroid(); 
  fOrgan["Thymus"] = new G4MIRDThymus(); 
  fOrgan["MaleGenitalia"] = new G4MIRDMaleGenitalia(); 
  fOrgan["Brain"] = new G4MIRDBrain(); 
  fOrgan["Stomach"] = new G4MIRDStomach();
  fOrgan["UpperLargeIntestine"] = new G4MIRDUpperLargeIntestine();
  fOrgan["LowerLargeIntestine"] = new G4MIRDLowerLargeIntestine();
  fOrgan["SmallIntestine"] = new G4MIRDSmallIntestine();
  fOrgan["Spleen"] = new G4MIRDSpleen(); 
  fOrgan["Pancreas"] = new G4MIRDPancreas();
  fOrgan["LeftKidney"] = new G4MIRDLeftKidney();
  fOrgan["RightKidney"] = new G4MIRDRightKidney();
  fOrgan["UrinaryBladder"] = new G4MIRDUrinaryBladder();
  fOrgan["Uterus"] = new G4MIRDUterus(); 
  fOrgan["Liver"] = new G4MIRDLiver(); 
  fOrgan["LeftLung"]= new G4MIRDLeftLung();
  fOrgan["RightLung"] = new G4MIRDRightLung();
  fOrgan["LeftOvary"] = new G4MIRDLeftOvary();
  fOrgan["RightOvary"] = new G4MIRDRightOvary();
  fOrgan["LeftTeste"] = new G4MIRDLeftTeste();
  fOrgan["RightTeste"] = new G4MIRDRightTeste();
  fOrgan["RightBreast"] = new G4MIRDRightBreast();
  fOrgan["LeftBreast"] = new G4MIRDLeftBreast();
  fOrgan["LeftAdrenal"]= new G4MIRDLeftAdrenal();  
  fOrgan["RightAdrenal"]= new G4MIRDRightAdrenal(); 
}

G4MIRDBodyFactory::~G4MIRDBodyFactory()
{
  delete fOrgan["Head"]; fOrgan["Head"]=nullptr;
  delete fOrgan["RightLeg"]; fOrgan["RightLeg"]=nullptr;
  delete fOrgan["LeftLeg"]; fOrgan["LeftLeg"]=nullptr;
  delete fOrgan["Trunk"]; fOrgan["Trunk"]=nullptr;
  delete fOrgan["RightScapula"];fOrgan["RightScapula"]=nullptr;
  delete fOrgan["LeftScapula"];fOrgan["LeftScapula"]=nullptr;
  delete fOrgan["RightLegBone"]; fOrgan["RightLegBone"]=nullptr;
  delete fOrgan["LeftLegBone"]; fOrgan["LeftLegBone"]=nullptr;
  delete fOrgan["RibCage"]; fOrgan["RibCage"] =nullptr;  
  delete fOrgan["MiddleLowerSpine"]; fOrgan["MidlleLowerSpine"]=nullptr;
  delete fOrgan["UpperSpine"]; fOrgan["UpperSpine"]=nullptr;
  delete fOrgan["Skull"]; fOrgan["Skull"]=nullptr;
  delete fOrgan["RightArmBone"]; fOrgan["RightArmBone"]=nullptr;
  delete fOrgan["LeftArmBone"]; fOrgan["LeftArmBone"]=nullptr;
  delete fOrgan["RightClavicle"]; fOrgan["RightClavicle"]=nullptr;
  delete fOrgan["LeftClavicle"]; fOrgan["LeftClavicle"]=nullptr;
  delete fOrgan["Pelvis"]; fOrgan["Pelvis"]=nullptr;  
  delete fOrgan["RightAdrenal"]; fOrgan["RightAdrenal"]=nullptr;
  delete fOrgan["LeftAdrenal"]; fOrgan["LeftAdrenal"]=nullptr;
  delete fOrgan["LeftBreast"]; fOrgan["LeftBreast"]=nullptr;
  delete fOrgan["RightBreast"]; fOrgan["RightBreast"]=nullptr;
  delete fOrgan["RightOvary"]; fOrgan["RightOvary"]=nullptr;
  delete fOrgan["LeftOvary"]; fOrgan["LeftOvary"]=nullptr;
  delete fOrgan["RightTeste"]; fOrgan["RightTeste"]=nullptr;
  delete fOrgan["LeftTeste"]; fOrgan["LeftTeste"]=nullptr;
  delete fOrgan["RightLung"]; fOrgan["RightLung"] =nullptr;
  delete fOrgan["LeftLung"]; fOrgan["LeftLung"]=nullptr;
  delete fOrgan["Uterus"]; fOrgan["Uterus"]=nullptr;
  delete fOrgan["UrinaryBladder"]; fOrgan["UrinaryBladder"]=nullptr;
  delete fOrgan["RightKidney"]; fOrgan["RightKidney"]=nullptr;  
  delete fOrgan["LeftKidney"]; fOrgan["LeftKidney"]=nullptr; 
  delete fOrgan["Pancreas"]; fOrgan["Pancreas"]=nullptr;  
  delete fOrgan["Spleen"]; fOrgan["Spleen"]=nullptr; 
  delete fOrgan["LowerLargeIntestine"]; fOrgan["LowerLargeIntestine"]=nullptr; 
  delete fOrgan["SmallIntestine"]; fOrgan["SmallIntestine"]=nullptr; 
  delete fOrgan["UpperLargeIntestine"]; fOrgan["UpperLargeIntestine"]=nullptr; 
  delete fOrgan["Stomach"]; fOrgan["Stomach"]=nullptr;  
  delete fOrgan["Brain"]; fOrgan["Brain"]=nullptr;
  delete fOrgan["Heart"]; fOrgan["Heart"]=nullptr;
  delete fOrgan["Thymus"]; fOrgan["Thymus"]=nullptr;
  delete fOrgan["MaleGenitalia"]; fOrgan["MaleGenitalia"]=nullptr;
  delete fOrgan["Thyroid"]; fOrgan["Thyroid"]=nullptr;
  delete fOrgan["Liver"]; fOrgan["Liver"]=nullptr;
}

G4VPhysicalVolume* G4MIRDBodyFactory::CreateOrgan(const G4String& organ_name,G4VPhysicalVolume* motherVolume,
						  const G4String& colourName, G4bool visAttribute,
						  G4bool sensitivity)
{
  return fOrgan[organ_name]->Construct(organ_name,motherVolume,colourName, visAttribute, sensitivity);
}


