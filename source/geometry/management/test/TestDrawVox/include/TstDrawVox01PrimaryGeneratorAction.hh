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

#ifndef TstDrawVox01PrimaryGeneratorAction_h
#define TstDrawVox01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4VisExtent.hh"

class G4ParticleGun;
class G4Event;
class TstDrawVox01PrimaryGeneratorMessenger;
class G4VPhysicalVolume;

class TstDrawVox01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  enum Action {
    standardGun,
    randomDirectionGun,
    randomPositionGun,
    randomPositionAndDirectionGun
  };

  TstDrawVox01PrimaryGeneratorAction();
  ~TstDrawVox01PrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);
  void SelectPrimaryGeneratorAction (Action);

private:
  Action generatorAction;
  G4ParticleGun* particleGun;
  TstDrawVox01PrimaryGeneratorMessenger* messenger;
  G4VPhysicalVolume* worldVolume;
  G4VisExtent worldExtent;
};

#endif
