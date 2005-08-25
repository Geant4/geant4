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
#ifndef G4VPhantomBuilder_h
#define G4VPhantomBuilder_h 1

#include "G4VPhysicalVolume.hh"
class G4VPhysicalVolume;
class G4VPhantomBuilder
{
public:

  G4VPhantomBuilder();
  virtual ~G4VPhantomBuilder();

  virtual void BuildWorld() = 0;

  virtual void BuildHead(G4bool) = 0;
  virtual void BuildTrunk(G4bool) = 0;
  virtual void BuildLegs(G4bool) = 0;
  virtual void BuildNeck(G4bool) = 0;

  virtual void BuildUpperSpine(G4bool) = 0;
  virtual void BuildMiddleLowerSpine(G4bool) = 0;
  virtual void BuildLegBone(G4bool) = 0;
  virtual void BuildArmBone(G4bool) = 0;
  virtual void BuildSkull(G4bool) = 0;
  virtual void BuildRibCage(G4bool) = 0;
  virtual void BuildPelvis(G4bool) = 0;
  virtual void BuildScapulae(G4bool) = 0;
  virtual void BuildClavicles(G4bool) = 0;

  virtual void BuildBrain(G4bool) = 0;

  virtual void BuildHeart(G4bool) = 0;

  virtual void BuildLung(G4bool) = 0;

  virtual void BuildStomach(G4bool) = 0;
  virtual void BuildUpperLargeIntestine(G4bool) = 0;
  virtual void BuildLowerLargeIntestine(G4bool) = 0;
  virtual void BuildEsophagus(G4bool) = 0;

  virtual void BuildKidney(G4bool) = 0;
  virtual void BuildAdrenal(G4bool) = 0;
  virtual void BuildLiver(G4bool) = 0;
  virtual void BuildPancreas(G4bool) = 0;
  virtual void BuildSpleen(G4bool) = 0;
  virtual void BuildUrinaryBladder(G4bool) = 0;

  virtual void BuildThyroid(G4bool) = 0;

  virtual void SetSex(G4String) = 0;
  virtual void SetModel(G4String) = 0;

};
#endif
