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
