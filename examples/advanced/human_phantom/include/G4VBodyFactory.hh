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

  virtual G4VPhysicalVolume* CreateHead(G4VPhysicalVolume*, G4String) = 0;
  virtual G4VPhysicalVolume* CreateTrunk(G4VPhysicalVolume*, G4String) = 0;
  virtual void CreateLegs(G4VPhysicalVolume*) = 0;
  virtual void CreateNeck(G4VPhysicalVolume*) = 0;

  virtual void CreateSpine(G4VPhysicalVolume*) = 0;
  virtual void CreateLegBone(G4VPhysicalVolume*) = 0;
  virtual void CreateArmBone(G4VPhysicalVolume*) = 0;
  virtual void CreateSkull(G4VPhysicalVolume*) = 0;
  virtual void CreateRibCage(G4VPhysicalVolume*) = 0;
  virtual void CreatePelvis(G4VPhysicalVolume*) = 0;
  virtual void CreateScapulae(G4VPhysicalVolume*) = 0;
  virtual void CreateClavicles(G4VPhysicalVolume*) = 0;

  virtual void CreateBreast(G4VPhysicalVolume*) = 0;
  virtual void CreateUterus(G4VPhysicalVolume*) = 0;
  virtual void CreateOvary(G4VPhysicalVolume*) = 0;

  virtual void CreateTestes(G4VPhysicalVolume*) = 0;
  virtual void CreateMaleGenitalia(G4VPhysicalVolume*) = 0;

  virtual void CreateBrain(G4VPhysicalVolume*) = 0;

  virtual void CreateHeart(G4VPhysicalVolume*) = 0;

  virtual void CreateLung(G4VPhysicalVolume*) = 0;

  virtual void CreateStomach(G4VPhysicalVolume*) = 0;
  virtual void CreateIntestine(G4VPhysicalVolume*) = 0;
  virtual void CreateEsophagus(G4VPhysicalVolume*) = 0;

  virtual void CreateKidney(G4VPhysicalVolume*) = 0;
  virtual void CreateAdrenal(G4VPhysicalVolume*) = 0;
  virtual void CreateLiver(G4VPhysicalVolume*) = 0;
  virtual void CreatePancreas(G4VPhysicalVolume*) = 0;
  virtual void CreateUrinaryBladder(G4VPhysicalVolume*) = 0;

  virtual void CreateThyroid(G4VPhysicalVolume*) = 0;


};
#endif
