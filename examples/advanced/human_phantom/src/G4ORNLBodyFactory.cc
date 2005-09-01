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
//
#include "G4ORNLBodyFactory.hh"
#include "G4ORNLStomach.hh"
#include "G4ORNLUpperLargeIntestine.hh"
#include "G4ORNLLowerLargeIntestine.hh"
#include "G4ORNLEsophagus.hh"
#include "G4ORNLKidney.hh"
#include "G4ORNLAdrenal.hh"
#include "G4ORNLLiver.hh"
#include "G4ORNLPancreas.hh"
#include "G4ORNLSpleen.hh"
#include "G4ORNLUrinaryBladder.hh"
#include "G4ORNLLung.hh"
#include "G4ORNLHeart.hh"
#include "G4ORNLBrain.hh"
#include "G4ORNLHead.hh"
#include "G4ORNLTrunk.hh"
#include "G4ORNLLegs.hh"
#include "G4ORNLThyroid.hh"
#include "G4ORNLUterus.hh"
#include "G4ORNLBreast.hh"
#include "G4ORNLOvary.hh"
#include "G4ORNLUpperSpine.hh"
#include "G4ORNLMiddleLowerSpine.hh"
#include "G4ORNLLegBone.hh"
#include "G4ORNLArmBone.hh"
#include "G4ORNLSkull.hh"
#include "G4ORNLRibCage.hh"
#include "G4ORNLPelvis.hh"
#include "G4ORNLTestes.hh"


G4ORNLBodyFactory::G4ORNLBodyFactory()
{
}

G4ORNLBodyFactory::~G4ORNLBodyFactory()
{
}

G4VPhysicalVolume* G4ORNLBodyFactory::CreateStomach(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLStomach* stomach = new G4ORNLStomach();
  return stomach -> ConstructStomach(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateUpperLargeIntestine(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLUpperLargeIntestine* upperLargeIntestine = new G4ORNLUpperLargeIntestine();
  return upperLargeIntestine -> ConstructUpperLargeIntestine(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateLowerLargeIntestine(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLLowerLargeIntestine* lowerLargeIntestine = new G4ORNLLowerLargeIntestine();
  return lowerLargeIntestine -> ConstructLowerLargeIntestine(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateEsophagus(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLEsophagus* esophagus = new G4ORNLEsophagus();
  return esophagus -> ConstructEsophagus(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateKidney(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLKidney* kidney = new G4ORNLKidney();
  return kidney -> ConstructKidney(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateAdrenal(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLAdrenal* adrenal = new G4ORNLAdrenal();
  return adrenal -> ConstructAdrenal(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateLiver(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLLiver* liver = new G4ORNLLiver();
  return liver -> ConstructLiver(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreatePancreas(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLPancreas* pancreas = new G4ORNLPancreas();
  return pancreas -> ConstructPancreas(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateSpleen(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLSpleen* spleen = new G4ORNLSpleen();
  return spleen -> ConstructSpleen(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateUrinaryBladder(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLUrinaryBladder* urinaryBladder = new G4ORNLUrinaryBladder();
  return urinaryBladder -> ConstructUrinaryBladder(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateLung(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLLung* lung = new G4ORNLLung();
  return lung -> ConstructLung(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateHeart(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLHeart* heart = new G4ORNLHeart();
  return heart -> ConstructHeart(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateBrain(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLBrain* brain = new G4ORNLBrain();
  return brain -> ConstructBrain(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateHead(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLHead* head = new G4ORNLHead();
  return head -> ConstructHead(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateTrunk(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLTrunk* trunk = new G4ORNLTrunk();
  return trunk -> ConstructTrunk(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateLegs(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLLegs* legs = new G4ORNLLegs();
  return legs -> ConstructLegs(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateThyroid(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLThyroid* thyroid = new G4ORNLThyroid();
  return thyroid -> ConstructThyroid(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateUterus(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLUterus* uterus = new G4ORNLUterus();
  return uterus -> ConstructUterus(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateBreast(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLBreast* breast = new G4ORNLBreast();
  return breast -> ConstructBreast(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateOvary(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLOvary* ovary = new G4ORNLOvary();
  return ovary -> ConstructOvary(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateUpperSpine(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLUpperSpine* upperSpine = new G4ORNLUpperSpine();
  return upperSpine -> ConstructUpperSpine(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateMiddleLowerSpine(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLMiddleLowerSpine* middleLowerSpine = new G4ORNLMiddleLowerSpine();
  return middleLowerSpine -> ConstructMiddleLowerSpine(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateLegBone(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLLegBone* legBone = new G4ORNLLegBone();
  return legBone -> ConstructLegBone(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateArmBone(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLArmBone* armBone = new G4ORNLArmBone();
  return armBone -> ConstructArmBone(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateSkull(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLSkull* skull = new G4ORNLSkull();
  return skull -> ConstructSkull(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateRibCage(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLRibCage* ribCage = new G4ORNLRibCage();
  return ribCage -> ConstructRibCage(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreatePelvis(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLPelvis* pelvis = new G4ORNLPelvis();
  return pelvis -> ConstructPelvis(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateTestes(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4ORNLTestes* testes = new G4ORNLTestes();
  return testes -> ConstructTestes(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateNeck(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
   G4cout << motherVolume << ">>>> Cannot create Neck of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl;
   return 0;
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateScapulae(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
   G4cout << motherVolume << ">>>> Cannot create Scapulae of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl;
   return 0;
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateClavicles(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
   G4cout << motherVolume << ">>>> Cannot create Clavicles of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl;
   return 0;
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateMaleGenitalia(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
   G4cout << motherVolume << ">>>> Cannot create Male Genitalia of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl;
   return 0;
}
