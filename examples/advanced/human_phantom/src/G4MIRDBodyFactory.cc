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
  // Map with name of the organ and pointer to the MIRDOrgan class
  //  organ["ParameterisedRightBreast"] = new G4ParameterisedRightBreast();
  //organ["ParameterisedLeftBreast"] = new G4ParameterisedLeftBreast();
  organ["Head"] = new G4MIRDHead();
  organ["Trunk"] = new G4MIRDTrunk(); 
  organ["LeftLeg"] = new G4MIRDLeftLeg();
  organ["RightLeg"] = new G4MIRDRightLeg();

  organ["Skull"] = new G4MIRDSkull();
  organ["LeftArmBone"] = new G4MIRDLeftArmBone();
  organ["RightArmBone"] = new G4MIRDRightArmBone();
  organ["UpperSpine"] = new G4MIRDUpperSpine();
  organ["MiddleLowerSpine"] = new G4MIRDMiddleLowerSpine();
  organ["Pelvis"]= new G4MIRDPelvis();
  organ["RibCage"] = new G4MIRDRibCage(); 
  organ["LeftClavicle"]= new G4MIRDLeftClavicle();
  organ["RightClavicle"] = new G4MIRDRightClavicle();
  organ["LeftLegBone"] = new G4MIRDLeftLegBone();
  organ["RightLegBone"] = new G4MIRDRightLegBone();
  organ["LeftScapula"]= new G4MIRDLeftScapula(); 
  organ["RightScapula"]= new G4MIRDRightScapula(); 

  organ["Heart"] = new G4MIRDHeart();
  organ["Thyroid"] = new G4MIRDThyroid(); 
  organ["Thymus"] = new G4MIRDThymus(); 
  organ["MaleGenitalia"] = new G4MIRDMaleGenitalia(); 
  organ["Brain"] = new G4MIRDBrain(); 
  organ["Stomach"] = new G4MIRDStomach();
  organ["UpperLargeIntestine"] = new G4MIRDUpperLargeIntestine();
  organ["LowerLargeIntestine"] = new G4MIRDLowerLargeIntestine();
  organ["SmallIntestine"] = new G4MIRDSmallIntestine();
  organ["Spleen"] = new G4MIRDSpleen(); 
  organ["Pancreas"] = new G4MIRDPancreas();
  organ["LeftKidney"] = new G4MIRDLeftKidney();
  organ["RightKidney"] = new G4MIRDRightKidney();
  organ["UrinaryBladder"] = new G4MIRDUrinaryBladder();
  organ["Uterus"] = new G4MIRDUterus(); 
  organ["Liver"] = new G4MIRDLiver(); 
  organ["LeftLung"]= new G4MIRDLeftLung();
  organ["RightLung"] = new G4MIRDRightLung();
  organ["LeftOvary"] = new G4MIRDLeftOvary();
  organ["RightOvary"] = new G4MIRDRightOvary();
  organ["LeftTeste"] = new G4MIRDLeftTeste();
  organ["RightTeste"] = new G4MIRDRightTeste();
  organ["RightBreast"] = new G4MIRDRightBreast();
  organ["LeftBreast"] = new G4MIRDLeftBreast();
  organ["LeftAdrenal"]= new G4MIRDLeftAdrenal();  
  organ["RightAdrenal"]= new G4MIRDRightAdrenal(); 
}

G4MIRDBodyFactory::~G4MIRDBodyFactory()
{
  delete organ["Head"]; organ["Head"]=0;
  delete organ["RightLeg"]; organ["RightLeg"]=0;
  delete organ["LeftLeg"]; organ["LeftLeg"]=0;
  delete organ["Trunk"]; organ["Trunk"]=0;

  delete organ["RightScapula"];organ["RightScapula"] =0;
  delete organ["LeftScapula"];organ["LeftScapula"] =0;
  delete organ["RightLegBone"]; organ["RightLegBone"]=0;
  delete organ["LeftLegBone"]; organ["LeftLegBone"]=0;
  delete organ["RibCage"]; organ["RibCage"] =0;  
  delete organ["MiddleLowerSpine"]; organ["MidlleLowerSpine"]=0;
  delete organ["UpperSpine"]; organ["UpperSpine"]=0;
  delete organ["Skull"]; organ["Skull"] =0;
  delete organ["RightArmBone"]; organ["RightArmBone"] =0;
  delete organ["LeftArmBone"]; organ["LeftArmBone"] =0;
  delete organ["RightClavicle"]; organ["RightClavicle"]=0;
  delete organ["LeftClavicle"]; organ["LeftClavicle"]=0;
  delete organ["Pelvis"]; organ["Pelvis"] =0;  

  delete organ["RightAdrenal"]; organ["RightAdrenal"]=0;
  delete organ["LeftAdrenal"]; organ["LeftAdrenal"]=0;
  delete organ["LeftBreast"]; organ["LeftBreast"]=0;
  delete organ["RightBreast"]; organ["RightBreast"]=0;
  delete organ["RightOvary"]; organ["RightOvary"]=0;
  delete organ["LeftOvary"]; organ["LeftOvary"]=0;
  delete organ["RightTeste"]; organ["RightTeste"]=0;
  delete organ["LeftTeste"]; organ["LeftTeste"]=0;
  delete organ["RightLung"]; organ["RightLung"] =0;
  delete organ["LeftLung"]; organ["LeftLung"]=0;
  delete organ["Uterus"]; organ["Uterus"]=0;
  delete organ["UrinaryBladder"]; organ["UrinaryBladder"]=0;
  delete organ["RightKidney"]; organ["RightKidney"] =0;  
  delete organ["LeftKidney"]; organ["LeftKidney"] =0; 
  delete organ["Pancreas"]; organ["Pancreas"] =0;  
  delete organ["Spleen"]; organ["Spleen"] =0; 
  delete organ["LowerLargeIntestine"]; organ["LowerLargeIntestine"] =0; 
  delete organ["SmallIntestine"]; organ["SmallIntestine"] =0; 
  delete organ["UpperLargeIntestine"]; organ["UpperLargeIntestine"] =0; 
  delete organ["Stomach"]; organ["Stomach"] =0;  
  delete organ["Brain"]; organ["Brain"]=0;
  delete organ["Heart"]; organ["Heart"]=0;
  delete organ["Thymus"]; organ["Thymus"]=0;
  delete organ["MaleGenitalia"]; organ["MaleGenitalia"]=0;
  delete organ["Thyroid"]; organ["Thyroid"]=0;
  delete organ["Liver"]; organ["Liver"]=0;
}



G4VPhysicalVolume* G4MIRDBodyFactory::CreateOrgan(const G4String& organ_name,G4VPhysicalVolume* motherVolume,
						  const G4String& colourName, G4bool visAttribute,
						  G4bool sensitivity)
{
  return organ[organ_name]->Construct(organ_name,motherVolume,colourName, visAttribute, sensitivity);
}


