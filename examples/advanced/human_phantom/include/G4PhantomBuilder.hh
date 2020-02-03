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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini INFN Perugia, Italy
//
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
  void BuildLeftScapula(const G4String&,G4bool,G4bool); 
  void BuildRightScapula(const G4String&,G4bool,G4bool);
  void BuildLeftAdrenal(const G4String&,G4bool,G4bool);  
  void BuildRightAdrenal(const G4String&,G4bool,G4bool); 
  void BuildLeftClavicle(const G4String&,G4bool,G4bool); 
  void BuildRightClavicle(const G4String&,G4bool,G4bool); 

  void BuildBrain(const G4String&,G4bool,G4bool);

  void BuildHeart(const G4String&,G4bool,G4bool);

  void BuildLeftLung(const G4String&,G4bool,G4bool);
  void BuildRightLung(const G4String&,G4bool,G4bool);

  void BuildStomach(const G4String&,G4bool,G4bool);
  void BuildSmallIntestine(const G4String&,G4bool,G4bool);
  void BuildUpperLargeIntestine(const G4String&,G4bool,G4bool);
  void BuildLowerLargeIntestine(const G4String&,G4bool,G4bool);

  void BuildLeftKidney(const G4String&,G4bool,G4bool);
  void BuildRightKidney(const G4String&,G4bool,G4bool);
  void BuildLiver(const G4String&,G4bool,G4bool);
  void BuildPancreas(const G4String&,G4bool,G4bool);
  void BuildSpleen(const G4String&,G4bool,G4bool);
  void BuildUrinaryBladder(const G4String&,G4bool,G4bool);

  void BuildThyroid(const G4String&,G4bool,G4bool);
  void BuildThymus(const G4String&,G4bool,G4bool);

  void SetModel(G4String);
  void SetMotherVolume(G4VPhysicalVolume*);

 
G4VPhysicalVolume* GetPhantom();

protected: 
  G4VBodyFactory* body;

  G4String model;

  G4VPhysicalVolume* motherVolume;
  G4VPhysicalVolume* headVolume;
  G4VPhysicalVolume* trunkVolume;
  G4VPhysicalVolume* leftLegVolume;
  G4VPhysicalVolume* rightLegVolume;  
  G4VPhysicalVolume* maleGenitaliaVolume;
};
#endif
