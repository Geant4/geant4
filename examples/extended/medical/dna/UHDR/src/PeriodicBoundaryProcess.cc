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
/// \file PeriodicBoundaryProcess.cc
/// \brief Implementation of the PeriodicBoundaryProcess class

/*
 * Based on 'g4pbc'.
 * Copyright (c) 2020 Amentum Pty Ltd
 * team@amentum.space
 * The original open-source version of this code
 * may be found at https://github.com/amentumspace/g4pbc
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
 * is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 */
#include "PeriodicBoundaryProcess.hh"

#include "LogicalVolumePeriodic.hh"
#include "ParticleChangeForPeriodic.hh"

#include "G4EventManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4Navigator.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"
#include "G4VTrajectory.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PeriodicBoundaryProcess::PeriodicBoundaryProcess(const G4String& processName, G4ProcessType type,
                                                 G4bool per_x, G4bool per_y, G4bool per_z)
  : G4VDiscreteProcess(processName, type), fPeriodicX(per_x), fPeriodicY(per_y), fPeriodicZ(per_z)
{
  fkCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  pParticleChange = &fParticleChange;
  G4PathFinder::GetInstance()->SetVerboseLevel(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* PeriodicBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  fTheStatus = Undefined;

  fParticleChange.InitializeForPostStep(aTrack);

  const G4Step* pStep = &aStep;

  G4bool isOnBoundary = (pStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);

  if (!isOnBoundary) {
    fTheStatus = NotAtBoundary;
    if (verboseLevel > 0) {
      BoundaryProcessVerbose();
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  if (aTrack.GetParentID() == 0)  // not primary particle (electron ?)
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  auto thePrePV = pStep->GetPreStepPoint()->GetPhysicalVolume();
  auto thePostPV = pStep->GetPostStepPoint()->GetPhysicalVolume();

  if (verboseLevel > 0) {
    G4cout << " Particle at Boundary! " << G4endl;
    if (thePrePV) {
      G4cout << " thePrePV:  " << thePrePV->GetName() << G4endl;
    }
    if (thePostPV) {
      G4cout << " thePostPV: " << thePostPV->GetName() << G4endl;
    }
    G4cout << "step length " << aTrack.GetStepLength() << G4endl;
    G4cout << "ParentID : " << aTrack.GetParentID() << "  TrackID : " << aTrack.GetTrackID()
           << " Position : " << aTrack.GetPosition()
           << "  aTrack Energy : " << G4BestUnit(aTrack.GetKineticEnergy(), "Energy") << G4endl;
  }

  // avoid trapped particles at boundaries by testing for minimum step length
  if (aTrack.GetStepLength() <= fkCarTolerance / 2) {
    fTheStatus = StepTooSmall;
    if (verboseLevel > 0) {
      BoundaryProcessVerbose();
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

  // store the current values
  fOldMomentum = aParticle->GetMomentumDirection();
  fOldPolarization = aParticle->GetPolarization();
  fOldPosition = pStep->GetPostStepPoint()->GetPosition();
  fNewPosition = fOldPosition;

  if (verboseLevel > 0) {
    G4cout << " Old Momentum Direction: " << fOldMomentum << G4endl;
    G4cout << " Old Position: " << fNewPosition << G4endl;
    G4cout << "Get aTrack.GetParentID() : " << aTrack.GetParentID() << G4endl;
  }

  auto theGlobalPoint = pStep->GetPostStepPoint()->GetPosition();

  G4bool valid = false;
  //  Use the new method for Exit Normal in global coordinates,
  //    which provides the normal more reliably.
  // get from the real-world navigator (process not applied to parallel worlds)
  fTheGlobalNormal = G4TransportationManager::GetTransportationManager()
                       ->GetNavigatorForTracking()
                       ->GetGlobalExitNormal(theGlobalPoint, &valid);

  if (valid) {
    fTheGlobalNormal = -fTheGlobalNormal;
  }
  else {
    G4cout << "global normal " << fTheGlobalNormal << G4endl;
    G4ExceptionDescription ed;
    ed << " PeriodicBoundaryProcess::PostStepDoIt(): "
       << " The Navigator reports that it returned an invalid normal" << G4endl;
    G4Exception("G4PeriodicBoundaryProcess::PostStepDoIt", "PerBoun01", FatalException, ed,
                "Invalid Surface Normal - Geometry must return valid surface normal");
  }

  G4bool isWrongDirection = fOldMomentum * fTheGlobalNormal > 0.0;
  if (isWrongDirection) {
    if (verboseLevel > 0) {
      G4cout << "ftheGlobalNormal points in a wrong direction." << G4endl;
      G4cout << "Invalid Surface Normal - Geometry must return valid surface \
        normal pointing in the right direction"
             << G4endl;
    }

    fTheGlobalNormal = -fTheGlobalNormal;
  }

  if (thePostPV == nullptr) {
    G4ExceptionDescription ed;
    ed << " PeriodicBoundaryProcess::PostStepDoIt(): "
       << " thePostPV == nullptr" << G4endl;
    G4Exception("G4PeriodicBoundaryProcess::PostStepDoIt", "PB14", FatalException, ed,
                "Invalid thePostPV");
  }
  auto lvol = thePostPV->GetLogicalVolume();

  if (lvol == nullptr) {
    G4ExceptionDescription ed;
    ed << " PeriodicBoundaryProcess::PostStepDoIt(): "
       << " lvol == nullptr" << G4endl;
    G4Exception("G4PeriodicBoundaryProcess::PostStepDoIt", "PB12", FatalException, ed,
                "Invalid lvol");
  }

  if (verboseLevel > 0) {
    G4cout << "Post step logical " << lvol->GetName() << G4endl;
  }

  G4LogicalVolume* dlvol = nullptr;

  if (lvol->GetNoDaughters() > 0) {
    if (verboseLevel > 0) {
      G4cout << "eldest daughter " << lvol->GetDaughter(0)->GetName() << G4endl;
    }

    dlvol = lvol->GetDaughter(0)->GetLogicalVolume();
  }

  if (dlvol && dlvol->IsExtended()) {
    if (verboseLevel > 0) {
      G4cout << " Logical surface, periodic " << G4endl;
    }

    G4bool on_x = fTheGlobalNormal.isParallel(G4ThreeVector(1, 0, 0));
    G4bool on_y = fTheGlobalNormal.isParallel(G4ThreeVector(0, 1, 0));
    G4bool on_z = fTheGlobalNormal.isParallel(G4ThreeVector(0, 0, 1));

    // make sure that we are at a plane
    G4bool on_plane = (on_x || on_y || on_z);

    if (!on_plane) {
      G4ExceptionDescription ed;
      ed << " G4PeriodicBoundaryProcess::ostStepDoIt(): "
         << " The particle is not on a surface of the cyclic world" << G4endl;
      G4Exception(
        "G4PeriodicBoundaryProcess::PostStepDoIt", "Periodic01", FatalException, ed,
        "Periodic boundary process must only occur for particle on periodic world surface");
    }
    else {
      G4bool on_x_and_periodic = (on_x && fPeriodicX);
      G4bool on_y_and_periodic = (on_y && fPeriodicY);
      G4bool on_z_and_periodic = (on_z && fPeriodicZ);

      G4bool on_a_periodic_plane = (on_x_and_periodic || on_y_and_periodic || on_z_and_periodic);

      if (on_a_periodic_plane) {
        if (verboseLevel > 0) {
          G4cout << " on periodic plane " << G4endl;
        }

        fTheStatus = Cycling;

        if (verboseLevel > 0) {
          G4cout << " periodic " << G4endl;
        }

        if (verboseLevel > 0) {
          G4cout << "Global normal " << fTheGlobalNormal << G4endl;
        }

        // translate a component of the position vector according to which plane we are on
        if (on_x_and_periodic) {
          fNewPosition.setX(-fNewPosition.x());
        }
        else if (on_y_and_periodic) {
          fNewPosition.setY(-fNewPosition.y());
        }
        else if (on_z_and_periodic) {
          fNewPosition.setZ(-fNewPosition.z());
        }
        else {
          G4cout << "global normal does not belong to periodic plane!!" << G4endl;
        }

        fNewMomentum = fOldMomentum.unit();
        fNewPolarization = fOldPolarization.unit();

        if (verboseLevel > 0) {
          G4cout << " New Position: " << fNewPosition << G4endl;
          G4cout << " New Momentum Direction: " << fNewMomentum << G4endl;
          G4cout << " New Polarization:       " << fNewPolarization << G4endl;
          BoundaryProcessVerbose();
        }

        fParticleChange.ProposeMomentumDirection(fNewMomentum);
        fParticleChange.ProposePolarization(fNewPolarization);
        fParticleChange.ProposePosition(fNewPosition);

        G4PathFinder::GetInstance()->ReLocate(fNewPosition);
        G4PathFinder::GetInstance()->ComputeSafety(fNewPosition);

        // we must notify the navigator that we have moved the particle artificially
        G4Navigator* gNavigator =
          G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

        // Inform the navigator that the previous Step calculated
        // by the geometry was taken in its entirety.
        gNavigator->SetGeometricallyLimitedStep();

        // Search the geometrical hierarchy for the volumes deepest in the hierarchy
        // containing the point in the global coordinate space.
        gNavigator->LocateGlobalPointAndSetup(fNewPosition, &fNewMomentum, false,
                                              false);  // do not ignore direction

        // Calculate the isotropic distance to the nearest boundary from the
        // specified point in the global coordinate system.
        gNavigator->ComputeSafety(fNewPosition);

        // Locates the volume containing the specified global point.
        // force drawing of the step prior to periodic the particle
        auto evtm = G4EventManager::GetEventManager();
        auto tckm = evtm->GetTrackingManager();
        auto pTrajectory = tckm->GimmeTrajectory();
        if (pTrajectory) {
          pTrajectory->AppendStep(pStep);
        }
      }
    }
  }
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PeriodicBoundaryProcess::GetMeanFreePath(const G4Track&, G4double,
                                                  G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PeriodicBoundaryProcess::BoundaryProcessVerbose()
{
  if (fStatusMessages.count(fTheStatus) > 0) {
    G4cout << PeriodicBoundaryProcess::fStatusMessages[fTheStatus] << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
