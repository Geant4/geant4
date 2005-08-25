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
// this is the abstract class which manages the body realisation
// in terms of geometry
//
#ifndef G4VBodyFactory_h
#define G4VBodyFactory_h 1

#include "G4VPhysicalVolume.hh"
class G4VPhysicalVolume;
class G4VBodyFactory
{
public:

  G4VBodyFactory();
  virtual ~G4VBodyFactory();

  virtual G4VPhysicalVolume* CreateHead(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateTrunk(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateLegs(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateNeck(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateUpperSpine(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateMiddleLowerSpine(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateLegBone(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateArmBone(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateSkull(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateRibCage(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreatePelvis(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateScapulae(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateClavicles(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateBreast(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateUterus(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateOvary(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateTestes(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateMaleGenitalia(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateBrain(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateHeart(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateLung(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateStomach(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateUpperLargeIntestine(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateLowerLargeIntestine(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateEsophagus(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateKidney(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateAdrenal(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateLiver(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreatePancreas(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateSpleen(G4VPhysicalVolume*,G4String,G4bool) = 0;
  virtual G4VPhysicalVolume* CreateUrinaryBladder(G4VPhysicalVolume*,G4String,G4bool) = 0;

  virtual G4VPhysicalVolume* CreateThyroid(G4VPhysicalVolume*,G4String,G4bool) = 0;


};
#endif
