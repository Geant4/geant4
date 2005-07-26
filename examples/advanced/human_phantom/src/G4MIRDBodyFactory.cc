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
// Code developed by:
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    GFactoryI.cc
//    *                             *
//    *******************************
//
// $Id: G4MIRDBodyFactory.cc,v 1.1 2005-07-26 13:41:22 gguerrie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4MIRDBodyFactory.hh"
#include "G4MIRDStomach.hh"
#include "G4MIRDIntestine.hh"
#include "G4MIRDKidney.hh"
#include "G4MIRDAdrenal.hh"
#include "G4MIRDLiver.hh"
#include "G4MIRDPancreas.hh"
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
#include "G4MIRDSpine.hh"
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

void G4MIRDBodyFactory::CreateStomach(G4VPhysicalVolume* motherVolume)
{
  G4MIRDStomach* stomach = new G4MIRDStomach();
  stomach -> ConstructStomach(motherVolume);
}

void G4MIRDBodyFactory::CreateIntestine(G4VPhysicalVolume* motherVolume)
{
  G4MIRDIntestine* intestine = new G4MIRDIntestine();
  intestine -> ConstructIntestine(motherVolume);
}

void G4MIRDBodyFactory::CreateKidney(G4VPhysicalVolume* motherVolume)
{
  G4MIRDKidney* kidney = new G4MIRDKidney();
  kidney -> ConstructKidney(motherVolume);
}
void G4MIRDBodyFactory::CreateAdrenal(G4VPhysicalVolume* motherVolume)
{
  G4MIRDAdrenal* adrenal = new G4MIRDAdrenal();
  adrenal -> ConstructAdrenal(motherVolume);
}
void G4MIRDBodyFactory::CreateLiver(G4VPhysicalVolume* motherVolume)
{
  G4MIRDLiver* liver = new G4MIRDLiver();
  liver -> ConstructLiver(motherVolume);
}
void G4MIRDBodyFactory::CreatePancreas(G4VPhysicalVolume* motherVolume)
{
  G4MIRDPancreas* pancreas = new G4MIRDPancreas();
  pancreas -> ConstructPancreas(motherVolume);
}
void G4MIRDBodyFactory::CreateUrinaryBladder(G4VPhysicalVolume* motherVolume)
{
  G4MIRDUrinaryBladder* urinaryBladder = new G4MIRDUrinaryBladder();
  urinaryBladder -> ConstructUrinaryBladder(motherVolume);
}
void G4MIRDBodyFactory::CreateLung(G4VPhysicalVolume* motherVolume)
{
  G4MIRDLung* lung = new G4MIRDLung();
  lung -> ConstructLung(motherVolume);
}
void G4MIRDBodyFactory::CreateHeart(G4VPhysicalVolume* motherVolume)
{
  G4MIRDHeart* heart = new G4MIRDHeart();
  heart -> ConstructHeart(motherVolume);
}
void G4MIRDBodyFactory::CreateBrain(G4VPhysicalVolume* motherVolume)
{
  G4MIRDBrain* brain = new G4MIRDBrain();
  brain -> ConstructBrain(motherVolume);
}

G4VPhysicalVolume* G4MIRDBodyFactory::CreateHead(G4VPhysicalVolume* motherVolume, G4String sex)
{
  G4MIRDHead* head = new G4MIRDHead();
  head -> ConstructHead(motherVolume, sex);
  return head->GetHead();
}
G4VPhysicalVolume* G4MIRDBodyFactory::CreateTrunk(G4VPhysicalVolume* motherVolume, G4String sex)
{
  G4MIRDTrunk* trunk = new G4MIRDTrunk();
  trunk -> ConstructTrunk(motherVolume,sex);
  return trunk -> GetTrunk();
}
void G4MIRDBodyFactory::CreateLegs(G4VPhysicalVolume* motherVolume)
{
  G4MIRDLegs* legs = new G4MIRDLegs();
  legs -> ConstructLegs(motherVolume);
}
void G4MIRDBodyFactory::CreateThyroid(G4VPhysicalVolume* motherVolume)
{
  G4MIRDThyroid* thyroid = new G4MIRDThyroid();
  thyroid -> ConstructThyroid(motherVolume);
}
void G4MIRDBodyFactory::CreateUterus(G4VPhysicalVolume* motherVolume)
{
  G4MIRDUterus* uterus = new G4MIRDUterus();
  uterus -> ConstructUterus(motherVolume);
}
void G4MIRDBodyFactory::CreateBreast(G4VPhysicalVolume* motherVolume)
{
  G4MIRDBreast* breast = new G4MIRDBreast();
  breast -> ConstructBreast(motherVolume);
}
void G4MIRDBodyFactory::CreateOvary(G4VPhysicalVolume* motherVolume)
{
  G4MIRDOvary* ovary = new G4MIRDOvary();
  ovary -> ConstructOvary(motherVolume);
}
void G4MIRDBodyFactory::CreateSpine(G4VPhysicalVolume* motherVolume)
{
  G4MIRDSpine* spine = new G4MIRDSpine();
  spine -> ConstructSpine(motherVolume);
}
void G4MIRDBodyFactory::CreateLegBone(G4VPhysicalVolume* motherVolume)
{
  G4MIRDLegBone* legBone = new G4MIRDLegBone();
  legBone -> ConstructLegBone(motherVolume);
}
void G4MIRDBodyFactory::CreateArmBone(G4VPhysicalVolume* motherVolume)
{
  G4MIRDArmBone* armBone = new G4MIRDArmBone();
  armBone -> ConstructArmBone(motherVolume);
}
void G4MIRDBodyFactory::CreateSkull(G4VPhysicalVolume* motherVolume)
{
  G4MIRDSkull* skull = new G4MIRDSkull();
  skull -> ConstructSkull(motherVolume);
}
void G4MIRDBodyFactory::CreateRibCage(G4VPhysicalVolume* motherVolume)
{
  G4MIRDRibCage* ribCage = new G4MIRDRibCage();
  ribCage -> ConstructRibCage(motherVolume);
}
void G4MIRDBodyFactory::CreatePelvis(G4VPhysicalVolume* motherVolume)
{
  G4MIRDPelvis* pelvis = new G4MIRDPelvis();
  pelvis -> ConstructPelvis(motherVolume);
}
void G4MIRDBodyFactory::CreateTestes(G4VPhysicalVolume* motherVolume)
{
  G4MIRDTestes* testes = new G4MIRDTestes();
  testes -> ConstructTestes(motherVolume);
}
void G4MIRDBodyFactory::CreateEsophagus(G4VPhysicalVolume* motherVolume)
{
}
void G4MIRDBodyFactory::CreateClavicles(G4VPhysicalVolume* motherVolume)
{
}
void G4MIRDBodyFactory::CreateScapulae(G4VPhysicalVolume* motherVolume)
{
}
void G4MIRDBodyFactory::CreateNeck(G4VPhysicalVolume* motherVolume)
{
}
void G4MIRDBodyFactory::CreateMaleGenitalia(G4VPhysicalVolume* motherVolume)
{
}
