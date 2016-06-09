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

#ifndef G4PhantomBuilder_h
#define G4PhantomBuilder_h 1

#include "G4VPhysicalVolume.hh"
#include "G4VPhantomBuilder.hh"
#include "globals.hh"

class G4VPhantomBuilder;
class G4VPhysicalVolume;
class G4VBodyFactory;
class G4PhantomBuilder: public G4VPhantomBuilder
{
public:
  G4PhantomBuilder();
  ~G4PhantomBuilder();

  void BuildWorld();

  void BuildHead(G4bool);
  void BuildTrunk(G4bool);
  void BuildLegs(G4bool);
  void BuildNeck(G4bool);

  void BuildUpperSpine(G4bool);
  void BuildMiddleLowerSpine(G4bool);
  void BuildLegBone(G4bool);
  void BuildArmBone(G4bool);
  void BuildSkull(G4bool);
//  void BuildRibCage(G4bool);
  void BuildPelvis(G4bool);
//  void BuildScapulae(G4bool);
//  void BuildClavicles(G4bool);

virtual void BuildBreast(G4bool)= 0;
virtual void BuildParameterisedBreast(G4bool) = 0;

virtual void BuildUterus(G4bool)=0;
virtual void BuildOvary(G4bool)=0;

virtual void BuildTestes(G4bool) =0;
virtual void BuildMaleGenitalia(G4bool)=0;

  void BuildBrain(G4bool);

  void BuildHeart(G4bool);

  void BuildLung(G4bool);

  void BuildStomach(G4bool);
  void BuildUpperLargeIntestine(G4bool);
  void BuildLowerLargeIntestine(G4bool);
 // void BuildEsophagus(G4bool);

  void BuildKidney(G4bool);
 // void BuildAdrenal(G4bool);
  void BuildLiver(G4bool);
  void BuildPancreas(G4bool);
  void BuildSpleen(G4bool);
  void BuildUrinaryBladder(G4bool);

  void BuildThyroid(G4bool);

  void SetSex(G4String);
  void SetModel(G4String);

  G4VPhysicalVolume* GetPhantom();

protected:
  G4String sex; 
  G4String model;
  G4bool sensitivity;
  G4VPhysicalVolume* motherVolume;
  G4VPhysicalVolume* headVolume;
  G4VPhysicalVolume* trunkVolume;
  G4VPhysicalVolume* legsVolume;
  G4VPhysicalVolume* neckVolume;
  G4VPhysicalVolume* maleGenitaliaVolume;
  G4VBodyFactory* body;
};
#endif
