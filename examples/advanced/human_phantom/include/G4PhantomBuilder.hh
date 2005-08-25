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
  void BuildRibCage(G4bool);
  void BuildPelvis(G4bool);
  void BuildScapulae(G4bool);
  void BuildClavicles(G4bool);

  void BuildBreast(G4bool);
  void BuildUterus(G4bool);
  void BuildOvary(G4bool);

  void BuildTestes(G4bool);
  void BuildMaleGenitalia(G4bool);

  void BuildBrain(G4bool);

  void BuildHeart(G4bool);

  void BuildLung(G4bool);

  void BuildStomach(G4bool);
  void BuildUpperLargeIntestine(G4bool);
  void BuildLowerLargeIntestine(G4bool);
  void BuildEsophagus(G4bool);

  void BuildKidney(G4bool);
  void BuildAdrenal(G4bool);
  void BuildLiver(G4bool);
  void BuildPancreas(G4bool);
  void BuildSpleen(G4bool);
  void BuildUrinaryBladder(G4bool);

  void BuildThyroid(G4bool);

  void SetSex(G4String);
  void SetModel(G4String);

  G4VPhysicalVolume* GetPhantom();

private:
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
