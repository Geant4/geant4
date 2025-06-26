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
// Author: Christian Velten (2025)

#include "G4MoleculeLocator.hh"

#include "G4IT.hh"
#include "G4ITTransportationManager.hh"
#include "G4Track.hh"
#include "G4TrackingInformation.hh"

G4ThreadLocal G4MoleculeLocator* G4MoleculeLocator::fpInstance = nullptr;

G4MoleculeLocator::G4MoleculeLocator()
{
  fNavigator = std::make_unique<G4ITNavigator>();
}

G4MoleculeLocator* G4MoleculeLocator::Instance()
{
  if (fpInstance == nullptr) {
    static G4ThreadLocalSingleton<G4MoleculeLocator> instance;
    fpInstance = instance.Instance();
  }
  if (!fpInstance->fIsInitialized) {
    fpInstance->Initialize();
  }
  return fpInstance;
}

void G4MoleculeLocator::Initialize()
{
  fNavigator->SetWorldVolume(G4ITTransportationManager::GetTransportationManager()
                               ->GetNavigatorForTracking()
                               ->GetWorldVolume());
  fIsInitialized = true;
}

G4TouchableHandle G4MoleculeLocator::LocateMoleculeTrack(const G4Track* pTrack)
{
  G4IT* pITrack = GetIT(pTrack);

  if (pITrack == nullptr) {
    G4Exception("G4MoleculeLocator::LocateMoleculeSetStateAndTouchable", "NOT_AN_IT", FatalErrorInArgument,
                "The track passed to this method appears to not hold an IT (molecule) object!");
  }

  std::unique_ptr<G4ITNavigatorState_Lock> tmpStateHolder;
  if (pITrack->GetTrackingInfo()->GetNavigatorState() != nullptr)
    fNavigator->SetNavigatorState(pITrack->GetTrackingInfo()->GetNavigatorState());
  else {
    fNavigator->NewNavigatorState();
    tmpStateHolder = std::unique_ptr<G4ITNavigatorState_Lock>(fNavigator->GetNavigatorState());
    // will be deleted once method goes out of scope
  }

  G4ThreeVector direction = pTrack->GetMomentumDirection();

  fNavigator->LocateGlobalPointAndSetup(pTrack->GetPosition(), &direction, false, false);

  G4TouchableHandle touchable = fNavigator->CreateTouchableHistory();
  return touchable;
}

void G4MoleculeLocator::LocateMoleculeSetStateAndTouchable(G4Track* pTrack)
{
  G4IT* pITrack = GetIT(pTrack);

  if (pITrack == nullptr) {
    G4Exception("G4MoleculeLocator::LocateMoleculeSetStateAndTouchable", "NOT_AN_IT", FatalErrorInArgument,
                "The track passed to this method appears to not hold an IT (molecule) object!");
  }

  fNavigator->NewNavigatorState();
  GetIT(pTrack)->GetTrackingInfo()->SetNavigatorState(fNavigator->GetNavigatorState());

  G4ThreeVector direction = pTrack->GetMomentumDirection();

  fNavigator->LocateGlobalPointAndSetup(pTrack->GetPosition(), &direction, false, false);

  G4TouchableHandle touchable = fNavigator->CreateTouchableHistory();
  pTrack->SetTouchableHandle(touchable);
  pTrack->SetNextTouchableHandle(touchable);
}
