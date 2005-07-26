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
  virtual void BuildHead() = 0;
  virtual void BuildTrunk() = 0;
//   virtual void BuildLegs(G4VPhysicalVolume*) = 0;
//   virtual void BuildNeck(G4VPhysicalVolume*) = 0;

//   virtual void BuildSpine(G4VPhysicalVolume*) = 0;
//   virtual void BuildLegBone(G4VPhysicalVolume*) = 0;
//   virtual void BuildArmBone(G4VPhysicalVolume*) = 0;
//   virtual void BuildSkull(G4VPhysicalVolume*) = 0;
//   virtual void BuildRibCage(G4VPhysicalVolume*) = 0;
//   virtual void BuildPelvis(G4VPhysicalVolume*) = 0;
//   virtual void BuildScapulae(G4VPhysicalVolume*) = 0;
//   virtual void BuildClavicles(G4VPhysicalVolume*) = 0;

//   virtual void BuildBreast(G4VPhysicalVolume*) = 0;
//   virtual void BuildUterus(G4VPhysicalVolume*) = 0;
//   virtual void BuildOvary(G4VPhysicalVolume*) = 0;

//   virtual void BuildTestes(G4VPhysicalVolume*) = 0;
//   virtual void BuildMaleGenitalia(G4VPhysicalVolume*) = 0;

//   virtual void BuildBrain(G4VPhysicalVolume*) = 0;

//   virtual void BuildHeart(G4VPhysicalVolume*) = 0;

//   virtual void BuildLung(G4VPhysicalVolume*) = 0;

//   virtual void BuildStomach(G4VPhysicalVolume*) = 0;
//   virtual void BuildIntestine(G4VPhysicalVolume*) = 0;
//   virtual void BuildEsophagus(G4VPhysicalVolume*) = 0;

//   virtual void BuildKidney(G4VPhysicalVolume*) = 0;
//   virtual void BuildAdrenal(G4VPhysicalVolume*) = 0;
//   virtual void BuildLiver(G4VPhysicalVolume*) = 0;
//   virtual void BuildPancreas(G4VPhysicalVolume*) = 0;
//   virtual void BuildUrinaryBladder(G4VPhysicalVolume*) = 0;

//   virtual void BuildThyroid(G4VPhysicalVolume*) = 0;
};
#endif
