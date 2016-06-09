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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
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
  virtual G4VPhysicalVolume* CreateParameterisedBreast(G4VPhysicalVolume*,G4String,G4bool) = 0;

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
