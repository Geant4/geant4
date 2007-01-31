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

#ifndef G4BasePhantomBuilder_h
#define G4BasePhantomBuilder_h 1

#include "G4VPhysicalVolume.hh"
class G4VPhysicalVolume;
class G4BasePhantomBuilder
{
public:

  G4BasePhantomBuilder();
  ~G4BasePhantomBuilder();
  // BuildHead(const G4String& volumeName, const G4String& colourName, G4bool visAttribute, G4bool) {return ;};
  void BuildHead(const G4String&,G4bool,G4bool) {return ;};
  void BuildTrunk(const G4String&,G4bool,G4bool) {return ;};
  void BuildLegs(const G4String&,G4bool,G4bool) {return ;};
  void BuildNeck(const G4String&,G4bool,G4bool) {return ;};

  void BuildUpperSpine(const G4String&,G4bool,G4bool) {return ;}
  void BuildMiddleLowerSpine(const G4String&,G4bool,G4bool) {return ;};
  void BuildLegBone(const G4String&,G4bool,G4bool) {return ;};
  void BuildLeftArmBone(const G4String&,G4bool,G4bool) {return ;}
  void BuildRightArmBone(const G4String&,G4bool,G4bool) {return ;}
  void BuildSkull(const G4String&,G4bool,G4bool) {return ;};
  void BuildRibCage(const G4String&,G4bool,G4bool) {return ;};
  void BuildPelvis(const G4String&,G4bool,G4bool) {return ;};
  //virtual void BuildScapulae(G4bool,G4bool) = 0;
  //virtual void BuildClavicles(G4bool,G4bool) = 0;

  void BuildBrain(const G4String&,G4bool,G4bool) {return ;};
  void BuildHeart(const G4String&,G4bool,G4bool) {return ;};
  void BuildLung(const G4String&,G4bool,G4bool) {return ;};
  void BuildStomach(const G4String&,G4bool,G4bool) {return ;};
  void BuildUpperLargeIntestine(const G4String&,G4bool,G4bool) {return ;};
  void BuildLowerLargeIntestine(const G4String&,G4bool,G4bool) {return ;};
  // virtual void BuildEsophagus(G4bool,G4bool) = 0;

  void BuildLeftKidney(const G4String&,G4bool,G4bool) {return ;};
  void BuildRightKidney(const G4String&,G4bool,G4bool) {return ;};
  // virtual void BuildAdrenal(G4bool,G4bool) = 0;
  void BuildLiver(const G4String&,G4bool,G4bool) {return ;};
  void BuildPancreas(const G4String&,G4bool,G4bool) {return ;};
  void BuildSpleen(const G4String&,G4bool,G4bool) {return ;};
  void BuildUrinaryBladder(const G4String& ,const G4String&,G4bool,G4bool) {return ;};
  void BuildThyroid(const G4String&,G4bool,G4bool) {return ;};

  void SetSex(G4String) {return ;};
  void SetModel(G4String) {return ;};
  void SetMotherVolume(G4VPhysicalVolume*) {return;};
  G4VPhysicalVolume* GetPhantom() {return 0;};

  virtual void BuildLeftOvary(const G4String&,G4bool,G4bool ) {return ;};
  virtual void BuildRightOvary(const G4String&,G4bool,G4bool) {return ;};
  virtual void BuildUterus(const G4String&,G4bool,G4bool){return;};
  virtual void BuildLeftBreast(const G4String&,G4bool,G4bool){return;};
  virtual void BuildRightBreast(const G4String&,G4bool,G4bool){return;};
};
#endif
