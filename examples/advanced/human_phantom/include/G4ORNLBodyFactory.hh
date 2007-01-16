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
#ifndef G4ORNLBodyFactory_h
#define G4ORNLBodyFactory_h 1

#include "G4VBodyFactory.hh"

class G4VBodyFactory;
class G4ORNLBodyFactory: public G4VBodyFactory
{
public:
  G4ORNLBodyFactory();
 ~G4ORNLBodyFactory();
  
  // G4VPhysicalVolume* CreateOrgan(G4VPhysicalVolume*,G4String,G4bool, G4String);

  G4VPhysicalVolume* CreateHead(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateTrunk(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateLegs(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateNeck(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateUpperSpine(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateMiddleLowerSpine(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateLegBone(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateLeftArmBone(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateRightArmBone(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateSkull(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateRibCage(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreatePelvis(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateScapulae(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateClavicles(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateBreast(G4VPhysicalVolume*,G4String,G4bool);
 G4VPhysicalVolume* CreateParameterisedBreast(G4VPhysicalVolume*,G4String,G4bool) { G4cout<< "Parameterised breast is not available for ORNL model"<< G4endl;
    return 0;};


  G4VPhysicalVolume* CreateUterus(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateOvary(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateTestes(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateMaleGenitalia(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateBrain(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateHeart(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateLung(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateStomach(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateUpperLargeIntestine(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateLowerLargeIntestine(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateEsophagus(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateLeftKidney(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateRightKidney(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateAdrenal(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateLiver(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreatePancreas(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateSpleen(G4VPhysicalVolume*,G4String,G4bool);
  G4VPhysicalVolume* CreateUrinaryBladder(G4VPhysicalVolume*,G4String,G4bool);

  G4VPhysicalVolume* CreateThyroid(G4VPhysicalVolume*,G4String,G4bool);


private:
  
};
#endif
