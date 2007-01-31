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

#ifndef G4PhantomBuilder_h
#define G4PhantomBuilder_h 1

#include "G4VPhysicalVolume.hh"
#include "G4BasePhantomBuilder.hh"
#include "globals.hh"

class G4BasePhantomBuilder;
class G4VPhysicalVolume;
class G4VBodyFactory;
class G4PhantomBuilder: public G4BasePhantomBuilder
{
public:
  G4PhantomBuilder();
  ~G4PhantomBuilder();

  void BuildHead(const G4String&,G4bool,G4bool);
  void BuildTrunk(const G4String&,G4bool,G4bool);
  void BuildLeftLeg(const G4String&,G4bool,G4bool);
  void BuildRightLeg(const G4String&,G4bool,G4bool);

  void BuildUpperSpine(const G4String&,G4bool,G4bool);
  void BuildMiddleLowerSpine(const G4String&,G4bool,G4bool);
  void BuildLeftLegBone(const G4String&,G4bool,G4bool);
  void BuildRightLegBone(const G4String&,G4bool,G4bool);
  void BuildLeftArmBone(const G4String&,G4bool,G4bool);
  void BuildRightArmBone(const G4String&,G4bool,G4bool);
  void BuildSkull(const G4String&,G4bool,G4bool);
  void BuildRibCage(const G4String&,G4bool,G4bool);
  void BuildPelvis(const G4String&,G4bool,G4bool);
//  void BuildScapulae(G4bool);
//  void BuildClavicles(G4bool);

  void BuildBrain(const G4String&,G4bool,G4bool);

  void BuildHeart(const G4String&,G4bool,G4bool);

  void BuildLeftLung(const G4String&,G4bool,G4bool);
  void BuildRightLung(const G4String&,G4bool,G4bool);

  void BuildStomach(const G4String&,G4bool,G4bool);
  void BuildUpperLargeIntestine(const G4String&,G4bool,G4bool);
  void BuildLowerLargeIntestine(const G4String&,G4bool,G4bool);
 // void BuildEsophagus(G4bool);

  void BuildLeftKidney(const G4String&,G4bool,G4bool);
  void BuildRightKidney(const G4String&,G4bool,G4bool);
 // void BuildAdrenal(G4bool);
  void BuildLiver(const G4String&,G4bool,G4bool);
  void BuildPancreas(const G4String&,G4bool,G4bool);
  void BuildSpleen(const G4String&,G4bool,G4bool);
  void BuildUrinaryBladder(const G4String&,G4bool,G4bool);

  void BuildThyroid(const G4String&,G4bool,G4bool);

  void SetModel(G4String);
  void SetMotherVolume(G4VPhysicalVolume*);

 
G4VPhysicalVolume* GetPhantom();

protected: 
  G4VBodyFactory* body;

  G4String model;
  G4bool sensitivity;
  G4VPhysicalVolume* motherVolume;
  G4VPhysicalVolume* headVolume;
  G4VPhysicalVolume* trunkVolume;
  G4VPhysicalVolume* leftLegVolume;
  G4VPhysicalVolume* rightLegVolume;  
  G4VPhysicalVolume* maleGenitaliaVolume;
};
#endif
