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

#include "G4FastSimHitMaker.hh"

#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHandle.hh"

#include "G4VFastSimSensitiveDetector.hh"

G4FastSimHitMaker::G4FastSimHitMaker()
{
  fTouchableHandle = new G4TouchableHistory();
  fpNavigator      = new G4Navigator();
  fNaviSetup       = false;
  fWorldWithSdName = "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FastSimHitMaker::~G4FastSimHitMaker() { delete fpNavigator; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4FastSimHitMaker::make(const G4FastHit& aHit, const G4FastTrack& aTrack)
{
  // do not make empty deposit
  if(aHit.GetEnergy() <= 0)
    return;
  // Locate the spot
  if(!fNaviSetup)
  {
    // Choose the world volume that contains the sensitive detector based on its
    // name (empty name for mass geometry)
    G4VPhysicalVolume* worldWithSD = nullptr;
    if(fWorldWithSdName.empty())
    {
      worldWithSD = G4TransportationManager::GetTransportationManager()
                      ->GetNavigatorForTracking()
                      ->GetWorldVolume();
    }
    else
    {
      worldWithSD =
        G4TransportationManager::GetTransportationManager()->GetParallelWorld(
          fWorldWithSdName);
    }
    fpNavigator->SetWorldVolume(worldWithSD);
    // use track global position
    fpNavigator->LocateGlobalPointAndUpdateTouchable(
      aTrack.GetPrimaryTrack()->GetPosition(), fTouchableHandle(), false);
    fNaviSetup = true;
  }
  else
  {
    // for further deposits use hit (local) position and local->global
    // transformation
    fpNavigator->LocateGlobalPointAndUpdateTouchable(
      aTrack.GetInverseAffineTransformation()->TransformPoint(
        aHit.GetPosition()),
      fTouchableHandle());
  }
  G4VPhysicalVolume* currentVolume = fTouchableHandle()->GetVolume();

  G4VSensitiveDetector* sensitive;
  if(currentVolume != 0)
  {
    sensitive = currentVolume->GetLogicalVolume()->GetSensitiveDetector();
    G4VFastSimSensitiveDetector* fastSimSensitive =
      dynamic_cast<G4VFastSimSensitiveDetector*>(sensitive);
    if(fastSimSensitive)
    {
      fastSimSensitive->Hit(&aHit, &aTrack, &fTouchableHandle);
    }
    else if(sensitive &&
            currentVolume->GetLogicalVolume()->GetFastSimulationManager())
    {
      G4cerr << "ERROR - G4FastSimHitMaker::make()" << G4endl
             << "        It is required to derive from the " << G4endl
             << "        G4VFastSimSensitiveDetector in " << G4endl
             << "        addition to the usual G4VSensitiveDetector class."
             << G4endl;
      G4Exception("G4FastSimHitMaker::make()", "InvalidSetup", FatalException,
                  "G4VFastSimSensitiveDetector interface not implemented.");
    }
  }
}
