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
  G4VPhysicalVolume* worldVolume;
  G4VPhysicalVolume* fDaughterVolume;
  G4VisExtent worldExtent;
  G4VisExtent fDaughterExtent;

};

///////////////////////////////////////////////////////////////////
//
//

inline void 
Tst01PrimaryGeneratorAction::SetGunPosition(G4ThreeVector pGun) 
{ 
  fGunPosition = pGun ;
 
  fPosition = sqrt( fGunPosition.x()*fGunPosition.x() +
                    fGunPosition.y()*fGunPosition.y() +
                    fGunPosition.z()*fGunPosition.z()     ) ;

  G4cout<<"fPosition = "<<fPosition<<G4endl ;

  // world extent

  fDaughterVolume = G4TransportationManager::GetTransportationManager ()
              -> GetNavigatorForTracking () -> GetWorldVolume ()
              ->GetLogicalVolume() -> GetDaughter(0)     ;
  if (fDaughterVolume) 
  {  
    fDaughterExtent = fDaughterVolume -> GetLogicalVolume () -> 
                                         GetSolid () -> GetExtent ();

    fSize = 0.6*sqrt( ( fDaughterExtent.GetXmax() - fDaughterExtent.GetXmin() )*
                  ( fDaughterExtent.GetXmax() - fDaughterExtent.GetXmin() ) +

                  ( fDaughterExtent.GetYmax() - fDaughterExtent.GetYmin() )*
                  ( fDaughterExtent.GetYmax() - fDaughterExtent.GetYmin() ) +

                  ( fDaughterExtent.GetZmax() - fDaughterExtent.GetZmin() )* 
                  ( fDaughterExtent.GetZmax() - fDaughterExtent.GetZmin() )     ) ;
  }
  else
  {
    fSize = 0.0 ;
  }
  G4cout<<"fSize = "<<fSize<<G4endl ;
}

#endif

//
//
/////////////////////////////////////////////////////////////////////
