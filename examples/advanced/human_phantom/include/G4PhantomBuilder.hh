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
  explicit G4PhantomBuilder();
  ~G4PhantomBuilder() override = default;

  void BuildHead(const G4String&,G4bool,G4bool) override;
  void BuildTrunk(const G4String&,G4bool,G4bool) override;
  void BuildLeftLeg(const G4String&,G4bool,G4bool) override;
  void BuildRightLeg(const G4String&,G4bool,G4bool) override;

  void BuildUpperSpine(const G4String&,G4bool,G4bool) override;
  void BuildMiddleLowerSpine(const G4String&,G4bool,G4bool) override;
  void BuildLeftLegBone(const G4String&,G4bool,G4bool) override;
  void BuildRightLegBone(const G4String&,G4bool,G4bool) override;
  void BuildLeftArmBone(const G4String&,G4bool,G4bool) override;
  void BuildRightArmBone(const G4String&,G4bool,G4bool) override;
  void BuildSkull(const G4String&,G4bool,G4bool) override;
  void BuildRibCage(const G4String&,G4bool,G4bool) override;
  void BuildPelvis(const G4String&,G4bool,G4bool) override;
  void BuildLeftScapula(const G4String&,G4bool,G4bool) override; 
  void BuildRightScapula(const G4String&,G4bool,G4bool) override;
  void BuildLeftAdrenal(const G4String&,G4bool,G4bool) override;  
  void BuildRightAdrenal(const G4String&,G4bool,G4bool) override; 
  void BuildLeftClavicle(const G4String&,G4bool,G4bool) override; 
  void BuildRightClavicle(const G4String&,G4bool,G4bool) override; 
  void BuildBrain(const G4String&,G4bool,G4bool) override;
  void BuildHeart(const G4String&,G4bool,G4bool) override;
  void BuildLeftLung(const G4String&,G4bool,G4bool) override;
  void BuildRightLung(const G4String&,G4bool,G4bool) override;
  void BuildStomach(const G4String&,G4bool,G4bool) override;
  void BuildSmallIntestine(const G4String&,G4bool,G4bool) override;
  void BuildUpperLargeIntestine(const G4String&,G4bool,G4bool) override;
  void BuildLowerLargeIntestine(const G4String&,G4bool,G4bool) override;
  void BuildLeftKidney(const G4String&,G4bool,G4bool) override;
  void BuildRightKidney(const G4String&,G4bool,G4bool) override;
  void BuildLiver(const G4String&,G4bool,G4bool) override;
  void BuildPancreas(const G4String&,G4bool,G4bool) override;
  void BuildSpleen(const G4String&,G4bool,G4bool) override;
  void BuildUrinaryBladder(const G4String&,G4bool,G4bool) override;
  void BuildThyroid(const G4String&,G4bool,G4bool) override;
  void BuildThymus(const G4String&,G4bool,G4bool) override;
  void SetModel(G4String) override;
  void SetMotherVolume(G4VPhysicalVolume*) override;

  G4VPhysicalVolume* GetPhantom() override;

protected: 
  G4VBodyFactory* fBody;
  G4String fModel;
  G4VPhysicalVolume* fMotherVolume;
  G4VPhysicalVolume* fHeadVolume;
  G4VPhysicalVolume* fTrunkVolume;
  G4VPhysicalVolume* fLeftLegVolume;
  G4VPhysicalVolume* fRightLegVolume;  
  G4VPhysicalVolume* fMaleGenitaliaVolume;
};
#endif
