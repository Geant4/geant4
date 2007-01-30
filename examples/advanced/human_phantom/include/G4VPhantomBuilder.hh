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

#ifndef G4VPhantomBuilder_h
#define G4VPhantomBuilder_h 1

#include "G4VPhysicalVolume.hh"
class G4VPhysicalVolume;
class G4VPhantomBuilder
{
public:

  G4VPhantomBuilder();
  ~G4VPhantomBuilder();

  void BuildHead(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildTrunk(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildLegs(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildNeck(G4bool, const G4String&, const G4String&, G4bool) {return ;};

  void BuildUpperSpine(G4bool, const G4String&, const G4String&, G4bool) {return ;}
  void BuildMiddleLowerSpine(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildLegBone(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildLeftArmBone(G4bool, const G4String&, const G4String&, G4bool) {return ;}
  void BuildRightArmBone(G4bool, const G4String&, const G4String&, G4bool) {return ;}
  void BuildSkull(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildRibCage(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildPelvis(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  //virtual void BuildScapulae(G4bool) = 0;
  //virtual void BuildClavicles(G4bool) = 0;

  void BuildBrain(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildHeart(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildLung(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildStomach(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildUpperLargeIntestine(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildLowerLargeIntestine(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  // virtual void BuildEsophagus(G4bool) = 0;

  void BuildLeftKidney(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildRightKidney(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  // virtual void BuildAdrenal(G4bool) = 0;
  void BuildLiver(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildPancreas(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildSpleen(G4bool,const G4String&, const G4String&, G4bool) {return ;};
  void BuildUrinaryBladder(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildThyroid(G4bool, const G4String&, const G4String&, G4bool) {return ;};
  void BuildLeftOvary(G4bool,const G4String&, const G4String&, G4bool ) {return ;};
  void BuildRightOvary(G4bool,const G4String&, const G4String&, G4bool) {return ;};
  void BuildUterus(G4bool,const G4String&, const G4String&, G4bool) {G4cout<< "BuildUterus!!!!!! in G4PhantomBuilder"<< G4endl;return ;};
  void BuildLeftBreast(G4bool,const G4String&, const G4String&, G4bool){return;};
  void BuildRightBreast(G4bool,const G4String&, const G4String&, G4bool){return;};
  //void BuildParameterisedBreast(G4bool) {return;};
  //void BuildTestes(G4bool) {return;};
  //void BuildMaleGenitalia(G4bool) {return;};

  void SetSex(G4String) {return ;};
  void SetModel(G4String) {return ;};
  void SetMotherVolume(G4VPhysicalVolume*) {return;};
};
#endif
