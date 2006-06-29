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

#ifndef Tst01PrimaryGeneratorAction_h
#define Tst01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VisExtent.hh"
#include "G4TransportationManager.hh"

class G4ParticleGun;
class G4Event;
class Tst01PrimaryGeneratorMessenger;
class G4VPhysicalVolume;

class Tst01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  enum Action { standardGun,
                randomDirectionGun,
                randomPositionGun,
                randomPositionAndDirectionGun,
                viewerGun,
                planeGun           } ;

  Tst01PrimaryGeneratorAction();
 ~Tst01PrimaryGeneratorAction();

public:

  void GeneratePrimaries(G4Event* anEvent);
  void SelectPrimaryGeneratorAction (Action);
  inline void SetGunPosition(G4ThreeVector pGun) ;

private:

  Action generatorAction;
  G4ParticleGun* particleGun;
  Tst01PrimaryGeneratorMessenger* messenger;
  G4ThreeVector fGunPosition ;
  G4double fPosition, fSize ;
};

///////////////////////////////////////////////////////////////////
//
//

inline void 
Tst01PrimaryGeneratorAction::SetGunPosition(G4ThreeVector pGun) 
{ 
  fGunPosition = pGun ;
  fPosition = std::sqrt( fGunPosition.x()*fGunPosition.x() +
                         fGunPosition.y()*fGunPosition.y() +
                         fGunPosition.z()*fGunPosition.z()  ) ;
  G4cout << "Absolute Gun Position = " << fPosition << G4endl;
}

#endif

//
//
/////////////////////////////////////////////////////////////////////
