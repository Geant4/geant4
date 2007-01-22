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
#include "G4MIRDBodyFactory.hh"
#include "G4MIRDStomach.hh"
#include "G4MIRDUpperLargeIntestine.hh"
#include "G4MIRDLowerLargeIntestine.hh"
#include "G4MIRDLeftKidney.hh"
#include "G4MIRDRightKidney.hh"
#include "G4MIRDAdrenal.hh"
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
#include "G4MIRDLegs.hh"
#include "G4MIRDThyroid.hh"
#include "G4MIRDUterus.hh"
#include "G4MIRDBreast.hh"
#include "G4MIRDRightOvary.hh"
#include "G4MIRDLeftOvary.hh"
#include "G4MIRDUpperSpine.hh"
#include "G4MIRDMiddleLowerSpine.hh"
#include "G4MIRDLegBone.hh"
#include "G4MIRDLeftArmBone.hh"
#include "G4MIRDRightArmBone.hh"
#include "G4MIRDSkull.hh"
#include "G4MIRDRibCage.hh"
#include "G4MIRDPelvis.hh"
#include "G4MIRDTestes.hh"
#include "G4ParameterisedBreast.hh"

G4MIRDBodyFactory::G4MIRDBodyFactory()
{
  // Map with name of the organ and pointer to the MIRDOrgan class
  organ["Head"] = new G4MIRDHead();
  organ["Trunk"] = new G4MIRDTrunk(); 
  organ["Legs"] = new G4MIRDLegs(); 
  organ["Brain"] = new G4MIRDBrain(); 
  organ["LeftArmBone"] = new G4MIRDLeftArmBone();
  organ["RightArmBone"] = new G4MIRDRightArmBone();
  organ["Skull"] = new G4MIRDSkull();
  organ["UpperSpine"] = new G4MIRDUpperSpine();
  organ["MiddleLowerSpine"] = new G4MIRDMiddleLowerSpine();
  organ["Pelvis"]= new G4MIRDPelvis();
  organ["Stomach"] = new G4MIRDStomach();
  organ["UpperLargeIntestine"] = new G4MIRDUpperLargeIntestine();
  organ["LowerLargeIntestine"] = new G4MIRDLowerLargeIntestine();
  organ["RibCage"] = new G4MIRDRibCage(); 
  organ["Spleen"] = new G4MIRDSpleen(); 
  organ["Pancreas"] = new G4MIRDPancreas();
  organ["LeftKidney"] = new G4MIRDLeftKidney();
  organ["RightKidney"] = new G4MIRDRightKidney();
  organ["UrinaryBladder"] = new G4MIRDUrinaryBladder();
  organ["Uterus"] = new G4MIRDUterus(); 
  organ["LeftLung"]= new G4MIRDLeftLung();
  organ["RightLung"] = new G4MIRDRightLung();
  organ["LeftOvary"] = new G4MIRDLeftOvary();
  organ["RightOvary"] = new G4MIRDRightOvary();
}

G4MIRDBodyFactory::~G4MIRDBodyFactory()
{ 
  delete organ["RightOvary"]; organ["RightOvary"]=0;
  delete organ["LeftOvary"]; organ["LeftOvary"]=0;
  delete organ["RightLung"]; organ["RightLung"] =0;
  delete organ["LeftLung"]; organ["LeftLung"]=0;
  delete organ["Uterus"]; organ["Uterus"]=0;
  delete organ["UrinaryBladder"]; organ["UrinaryBladder"]=0;
  delete organ["RightKidney"]; organ["RightKidney"] =0;  
  delete organ["LeftKidney"]; organ["LeftKidney"] =0; 
  delete organ["Pancreas"]; organ["Pancreas"] =0;  
  delete organ["Spleen"]; organ["Spleen"] =0; 
  delete organ["RibCage"]; organ["RibCage"] =0;  
  delete organ["LowerLargeIntestine"]; organ["LowerLargeIntestine"] =0; 
  delete organ["UpperLargeIntestine"]; organ["UpperLargeIntestine"] =0; 
  delete organ["Stomach"]; organ["Stomach"] =0;  
  delete organ["Pelvis"]; organ["Pelvis"] =0;  
  delete organ["MiddleLowerSpine"]; organ["MidlleLowerSpine"]=0;
  delete organ["UpperSpine"]; organ["UpperSpine"]=0;
  delete organ["Skull"]; organ["Skull"] =0;
  delete organ["RightArmBone"]; organ["RightArmBone"] =0;
  delete organ["LeftArmBone"]; organ["LeftArmBone"] =0;
  delete organ["Brain"]; organ["Brain"]=0;
  delete organ["Head"]; organ["Head"]=0;
  delete organ["Trunk"]; organ["Trunk"]=0;
  delete organ["Legs"]; organ["Legs"]=0;
}

G4VPhysicalVolume* G4MIRDBodyFactory::CreateOrgan(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity, G4String organ_name,
						  G4String logicalVolumeName, G4String colourName, G4bool solidColour)
{
 return organ[organ_name] -> ConstructOrgan(motherVolume, sex, sensitivity, organ_name, logicalVolumeName, colourName, solidColour);
}


