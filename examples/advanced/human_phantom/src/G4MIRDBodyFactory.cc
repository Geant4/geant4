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
#include "G4MIRDBodyFactory.hh"
#include "G4MIRDStomach.hh"
#include "G4MIRDUpperLargeIntestine.hh"
#include "G4MIRDLowerLargeIntestine.hh"
#include "G4MIRDKidney.hh"
#include "G4MIRDAdrenal.hh"
#include "G4MIRDLiver.hh"
#include "G4MIRDPancreas.hh"
#include "G4MIRDSpleen.hh"
#include "G4MIRDUrinaryBladder.hh"
#include "G4MIRDLung.hh"
#include "G4MIRDHeart.hh"
#include "G4MIRDBrain.hh"
#include "G4MIRDHead.hh"
#include "G4MIRDTrunk.hh"
#include "G4MIRDLegs.hh"
#include "G4MIRDThyroid.hh"
#include "G4MIRDUterus.hh"
#include "G4MIRDBreast.hh"
#include "G4MIRDOvary.hh"
#include "G4MIRDUpperSpine.hh"
#include "G4MIRDMiddleLowerSpine.hh"
#include "G4MIRDLegBone.hh"
#include "G4MIRDArmBone.hh"
#include "G4MIRDSkull.hh"
#include "G4MIRDRibCage.hh"
#include "G4MIRDPelvis.hh"
#include "G4MIRDTestes.hh"


G4MIRDBodyFactory::G4MIRDBodyFactory()
{
}

G4MIRDBodyFactory::~G4MIRDBodyFactory()
{
}

G4VPhysicalVolume* G4MIRDBodyFactory::CreateStomach(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDStomach* stomach = new G4MIRDStomach();
  return stomach -> ConstructStomach(motherVolume,sex,sensitivity);
}

G4VPhysicalVolume* G4MIRDBodyFactory::CreateUpperLargeIntestine(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDUpperLargeIntestine* upperLargeIntestine = new G4MIRDUpperLargeIntestine();
  return upperLargeIntestine -> ConstructUpperLargeIntestine(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateLowerLargeIntestine(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDLowerLargeIntestine* lowerLargeIntestine = new G4MIRDLowerLargeIntestine();
  return lowerLargeIntestine -> ConstructLowerLargeIntestine(motherVolume,sex,sensitivity);
}

G4VPhysicalVolume* G4MIRDBodyFactory::CreateKidney(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDKidney* kidney = new G4MIRDKidney();
  return kidney -> ConstructKidney(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateAdrenal(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDAdrenal* adrenal = new G4MIRDAdrenal();
  return adrenal -> ConstructAdrenal(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateLiver(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDLiver* liver = new G4MIRDLiver();
  return liver -> ConstructLiver(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreatePancreas(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDPancreas* pancreas = new G4MIRDPancreas();
  return pancreas -> ConstructPancreas(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateSpleen(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDSpleen* spleen = new G4MIRDSpleen();
  return spleen -> ConstructSpleen(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateUrinaryBladder(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDUrinaryBladder* urinaryBladder = new G4MIRDUrinaryBladder();
  return urinaryBladder -> ConstructUrinaryBladder(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateLung(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDLung* lung = new G4MIRDLung();
  return lung -> ConstructLung(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateHeart(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDHeart* heart = new G4MIRDHeart();
  return heart -> ConstructHeart(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateBrain(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDBrain* brain = new G4MIRDBrain();
  return brain -> ConstructBrain(motherVolume,sex,sensitivity);
}

G4VPhysicalVolume* G4MIRDBodyFactory::CreateHead(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDHead* head = new G4MIRDHead();
  return head -> ConstructHead(motherVolume, sex, sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateTrunk(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDTrunk* trunk = new G4MIRDTrunk();
  return trunk -> ConstructTrunk(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateLegs(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDLegs* legs = new G4MIRDLegs();
  return legs -> ConstructLegs(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateThyroid(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDThyroid* thyroid = new G4MIRDThyroid();
  return thyroid -> ConstructThyroid(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateUterus(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDUterus* uterus = new G4MIRDUterus();
  return uterus -> ConstructUterus(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateBreast(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDBreast* breast = new G4MIRDBreast();
  return breast -> ConstructBreast(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateOvary(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDOvary* ovary = new G4MIRDOvary();
  return ovary -> ConstructOvary(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateUpperSpine(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDUpperSpine* upperSpine = new G4MIRDUpperSpine();
  return upperSpine -> ConstructUpperSpine(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateMiddleLowerSpine(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDMiddleLowerSpine* middleLowerSpine = new G4MIRDMiddleLowerSpine();
  return middleLowerSpine -> ConstructMiddleLowerSpine(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateLegBone(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDLegBone* legBone = new G4MIRDLegBone();
  return legBone -> ConstructLegBone(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateArmBone(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDArmBone* armBone = new G4MIRDArmBone();
  return armBone -> ConstructArmBone(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateSkull(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDSkull* skull = new G4MIRDSkull();
  return skull -> ConstructSkull(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateRibCage(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDRibCage* ribCage = new G4MIRDRibCage();
  return ribCage -> ConstructRibCage(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreatePelvis(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDPelvis* pelvis = new G4MIRDPelvis();
  return pelvis -> ConstructPelvis(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateTestes(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4MIRDTestes* testes = new G4MIRDTestes();
  return testes -> ConstructTestes(motherVolume,sex,sensitivity);
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateEsophagus(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
  G4cout << motherVolume << ">>>> Cannot create Esophagus of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl; 
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateClavicles(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
   G4cout << motherVolume << ">>>> Cannot create Clavicles of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl; 
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateScapulae(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
   G4cout << motherVolume << ">>>> Cannot create Scapulae of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl; 
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateNeck(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
   G4cout << motherVolume << ">>>> Cannot create Neck of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl; 
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateMaleGenitalia(G4VPhysicalVolume* motherVolume, G4String sex, G4bool sensitivity)
{
   G4cout << motherVolume << ">>>> Cannot create Male Genitalia of " << sex 
	 << "phantom and with sensitivity = " << sensitivity 
	 << ": not yet implemented!" << G4endl; 
}
