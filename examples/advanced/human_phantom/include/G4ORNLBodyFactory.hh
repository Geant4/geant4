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
#ifndef G4ORNLBodyFactory_h
#define G4ORNLBodyFactory_h 1

#include "G4VBodyFactory.hh"

class G4VBodyFactory;
class G4ORNLBodyFactory: public G4VBodyFactory
{
public:
  G4ORNLBodyFactory();
 ~G4ORNLBodyFactory();

  G4VPhysicalVolume* CreateHead(G4VPhysicalVolume*, G4String);
  G4VPhysicalVolume* CreateTrunk(G4VPhysicalVolume*, G4String);
  void CreateLegs(G4VPhysicalVolume*);
  void CreateNeck(G4VPhysicalVolume*);

  void CreateSpine(G4VPhysicalVolume*);
  void CreateLegBone(G4VPhysicalVolume*);
  void CreateArmBone(G4VPhysicalVolume*);
  void CreateSkull(G4VPhysicalVolume*);
  void CreateRibCage(G4VPhysicalVolume*);
  void CreatePelvis(G4VPhysicalVolume*);
  void CreateScapulae(G4VPhysicalVolume*);
  void CreateClavicles(G4VPhysicalVolume*);

  void CreateBreast(G4VPhysicalVolume*);
  void CreateUterus(G4VPhysicalVolume*);
  void CreateOvary(G4VPhysicalVolume*);

  void CreateTestes(G4VPhysicalVolume*);
  void CreateMaleGenitalia(G4VPhysicalVolume*);

  void CreateBrain(G4VPhysicalVolume*);

  void CreateHeart(G4VPhysicalVolume*);

  void CreateLung(G4VPhysicalVolume*);

  void CreateStomach(G4VPhysicalVolume*);
  void CreateIntestine(G4VPhysicalVolume*);
  void CreateEsophagus(G4VPhysicalVolume*);

  void CreateKidney(G4VPhysicalVolume*);
  void CreateAdrenal(G4VPhysicalVolume*);
  void CreateLiver(G4VPhysicalVolume*);
  void CreatePancreas(G4VPhysicalVolume*);
  void CreateUrinaryBladder(G4VPhysicalVolume*);

  void CreateThyroid(G4VPhysicalVolume*);


private:
  
};
#endif
