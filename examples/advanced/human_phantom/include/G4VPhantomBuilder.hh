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
  //virtual void BuildRibCage(G4bool) = 0;
  virtual void BuildPelvis(G4bool) = 0;
  //virtual void BuildScapulae(G4bool) = 0;
  //virtual void BuildClavicles(G4bool) = 0;

  virtual void BuildBrain(G4bool) = 0;

  virtual void BuildHeart(G4bool) = 0;

  virtual void BuildLung(G4bool) = 0;

  virtual void BuildStomach(G4bool) = 0;
  virtual void BuildUpperLargeIntestine(G4bool) = 0;
  virtual void BuildLowerLargeIntestine(G4bool) = 0;
  //virtual void BuildEsophagus(G4bool) = 0;

  virtual void BuildKidney(G4bool) = 0;
//  virtual void BuildAdrenal(G4bool) = 0;
  virtual void BuildLiver(G4bool) = 0;
  virtual void BuildPancreas(G4bool) = 0;
  virtual void BuildSpleen(G4bool) = 0;
  virtual void BuildUrinaryBladder(G4bool) = 0;

  virtual void BuildThyroid(G4bool) = 0;

  virtual void SetSex(G4String) = 0;
  virtual void SetModel(G4String) = 0;

};
#endif
