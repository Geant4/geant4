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
// Previous authors: G. Guerrieri, S. Guatelli, and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//

#ifndef G4BasePhantomBuilder_h
#define G4BasePhantomBuilder_h 1

#include "G4VPhysicalVolume.hh"
class G4VPhysicalVolume;
class G4BasePhantomBuilder
{
public:

  G4BasePhantomBuilder();
  virtual ~G4BasePhantomBuilder();
 
  virtual void BuildHead(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildTrunk(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLeftLeg(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildRightLeg(const G4String&,G4bool,G4bool) {return ;};

  virtual void BuildUpperSpine(const G4String&,G4bool,G4bool) {return ;}
  virtual void BuildMiddleLowerSpine(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLeftLegBone(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildRightLegBone(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLeftArmBone(const G4String&,G4bool,G4bool) {return ;}
  virtual void BuildRightArmBone(const G4String&,G4bool,G4bool) {return ;}
  virtual void BuildSkull(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildRibCage(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildPelvis(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLeftScapula(const G4String&,G4bool,G4bool){return;};
  virtual void BuildRightScapula(const G4String&,G4bool,G4bool){return;};
  virtual void BuildLeftClavicle(const G4String&,G4bool,G4bool){return;};
  virtual void BuildRightClavicle(const G4String&,G4bool,G4bool){return;};
 
  virtual void BuildBrain(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildHeart(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLeftLung(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildRightLung(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildStomach(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildSmallIntestine(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildUpperLargeIntestine(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLowerLargeIntestine(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLeftKidney(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildRightKidney(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLeftAdrenal(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildRightAdrenal(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLiver(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildPancreas(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildSpleen(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildUrinaryBladder(const G4String& ,G4bool,G4bool) {return ;};
  virtual void BuildThyroid(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildThymus(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildLeftOvary(const G4String&,G4bool,G4bool ) {return ;};
  virtual void BuildRightOvary(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildUterus(const G4String&,G4bool,G4bool){return;};
  virtual void BuildLeftBreast(const G4String&,G4bool,G4bool){return;};
  virtual void BuildRightBreast(const G4String&,G4bool,G4bool){return;};
  virtual void BuildVoxelLeftBreast(const G4String&,G4bool,G4bool){return;};
  virtual void BuildVoxelRightBreast(const G4String&,G4bool,G4bool){return;};
  virtual void BuildMaleGenitalia(const G4String&,G4bool,G4bool){return;};
  virtual void BuildLeftTeste(const G4String&,G4bool,G4bool){return;};
  virtual void BuildRightTeste(const G4String&,G4bool,G4bool){return;};

  virtual void SetModel(G4String) {return ;};
  virtual void SetMotherVolume(G4VPhysicalVolume*) {return;};
  virtual G4VPhysicalVolume* GetPhantom() {return 0;};

};
#endif
