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
#include "G4ORNLIntestine.hh"
#include "G4ORNLEsophagus.hh"
#include "G4ORNLKidney.hh"
#include "G4ORNLAdrenal.hh"
#include "G4ORNLLiver.hh"
#include "G4ORNLPancreas.hh"
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
#include "G4ORNLSpine.hh"
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

void G4ORNLBodyFactory::CreateStomach(G4VPhysicalVolume* motherVolume)
{
  G4ORNLStomach* stomach = new G4ORNLStomach();
  stomach -> ConstructStomach(motherVolume);
}

void G4ORNLBodyFactory::CreateIntestine(G4VPhysicalVolume* motherVolume)
{
  G4ORNLIntestine* intestine = new G4ORNLIntestine();
  intestine -> ConstructIntestine(motherVolume);
}

void G4ORNLBodyFactory::CreateEsophagus(G4VPhysicalVolume* motherVolume)
{
  G4ORNLEsophagus* esophagus = new G4ORNLEsophagus();
  esophagus -> ConstructEsophagus(motherVolume);
}

void G4ORNLBodyFactory::CreateKidney(G4VPhysicalVolume* motherVolume)
{
  G4ORNLKidney* kidney = new G4ORNLKidney();
  kidney -> ConstructKidney(motherVolume);
}
void G4ORNLBodyFactory::CreateAdrenal(G4VPhysicalVolume* motherVolume)
{
  G4ORNLAdrenal* adrenal = new G4ORNLAdrenal();
  adrenal -> ConstructAdrenal(motherVolume);
}
void G4ORNLBodyFactory::CreateLiver(G4VPhysicalVolume* motherVolume)
{
  G4ORNLLiver* liver = new G4ORNLLiver();
  liver -> ConstructLiver(motherVolume);
}
void G4ORNLBodyFactory::CreatePancreas(G4VPhysicalVolume* motherVolume)
{
  G4ORNLPancreas* pancreas = new G4ORNLPancreas();
  pancreas -> ConstructPancreas(motherVolume);
}
void G4ORNLBodyFactory::CreateUrinaryBladder(G4VPhysicalVolume* motherVolume)
{
  G4ORNLUrinaryBladder* urinaryBladder = new G4ORNLUrinaryBladder();
  urinaryBladder -> ConstructUrinaryBladder(motherVolume);
}
void G4ORNLBodyFactory::CreateLung(G4VPhysicalVolume* motherVolume)
{
  G4ORNLLung* lung = new G4ORNLLung();
  lung -> ConstructLung(motherVolume);
}
void G4ORNLBodyFactory::CreateHeart(G4VPhysicalVolume* motherVolume)
{
  G4ORNLHeart* heart = new G4ORNLHeart();
  heart -> ConstructHeart(motherVolume);
}
void G4ORNLBodyFactory::CreateBrain(G4VPhysicalVolume* motherVolume)
{
  G4ORNLBrain* brain = new G4ORNLBrain();
  brain -> ConstructBrain(motherVolume);
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateHead(G4VPhysicalVolume* motherVolume, G4String sex)
{
  G4ORNLHead* head = new G4ORNLHead();
  head -> ConstructHead(motherVolume, sex);
  return head -> GetHead();
}
G4VPhysicalVolume* G4ORNLBodyFactory::CreateTrunk(G4VPhysicalVolume* motherVolume, G4String sex)
{
  G4ORNLTrunk* trunk = new G4ORNLTrunk();
  trunk -> ConstructTrunk(motherVolume, sex);
  return trunk -> GetTrunk(); 
}
void G4ORNLBodyFactory::CreateLegs(G4VPhysicalVolume* motherVolume)
{
  G4ORNLLegs* legs = new G4ORNLLegs();
  legs -> ConstructLegs(motherVolume);
}
void G4ORNLBodyFactory::CreateThyroid(G4VPhysicalVolume* motherVolume)
{
  G4ORNLThyroid* thyroid = new G4ORNLThyroid();
  thyroid -> ConstructThyroid(motherVolume);
}
void G4ORNLBodyFactory::CreateUterus(G4VPhysicalVolume* motherVolume)
{
  G4ORNLUterus* uterus = new G4ORNLUterus();
  uterus -> ConstructUterus(motherVolume);
}
void G4ORNLBodyFactory::CreateBreast(G4VPhysicalVolume* motherVolume)
{
  G4ORNLBreast* breast = new G4ORNLBreast();
  breast -> ConstructBreast(motherVolume);
}
void G4ORNLBodyFactory::CreateOvary(G4VPhysicalVolume* motherVolume)
{
  G4ORNLOvary* ovary = new G4ORNLOvary();
  ovary -> ConstructOvary(motherVolume);
}
void G4ORNLBodyFactory::CreateSpine(G4VPhysicalVolume* motherVolume)
{
  G4ORNLSpine* spine = new G4ORNLSpine();
  spine -> ConstructSpine(motherVolume);
}
void G4ORNLBodyFactory::CreateLegBone(G4VPhysicalVolume* motherVolume)
{
  G4ORNLLegBone* legBone = new G4ORNLLegBone();
  legBone -> ConstructLegBone(motherVolume);
}
void G4ORNLBodyFactory::CreateArmBone(G4VPhysicalVolume* motherVolume)
{
  G4ORNLArmBone* armBone = new G4ORNLArmBone();
  armBone -> ConstructArmBone(motherVolume);
}
void G4ORNLBodyFactory::CreateSkull(G4VPhysicalVolume* motherVolume)
{
  G4ORNLSkull* skull = new G4ORNLSkull();
  skull -> ConstructSkull(motherVolume);
}
void G4ORNLBodyFactory::CreateRibCage(G4VPhysicalVolume* motherVolume)
{
  G4ORNLRibCage* ribCage = new G4ORNLRibCage();
  ribCage -> ConstructRibCage(motherVolume);
}
void G4ORNLBodyFactory::CreatePelvis(G4VPhysicalVolume* motherVolume)
{
  G4ORNLPelvis* pelvis = new G4ORNLPelvis();
  pelvis -> ConstructPelvis(motherVolume);
}
void G4ORNLBodyFactory::CreateTestes(G4VPhysicalVolume* motherVolume)
{
  G4ORNLTestes* testes = new G4ORNLTestes();
  testes -> ConstructTestes(motherVolume);
}
void G4ORNLBodyFactory::CreateNeck(G4VPhysicalVolume* motherVolume)
{
}
void G4ORNLBodyFactory::CreateScapulae(G4VPhysicalVolume* motherVolume)
{
}
void G4ORNLBodyFactory::CreateClavicles(G4VPhysicalVolume* motherVolume)
{
}
void G4ORNLBodyFactory::CreateMaleGenitalia(G4VPhysicalVolume* motherVolume)
{
}
