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
//
/// \file optical/wls/src/WLSSteppingAction.cc
/// \brief Implementation of the WLSSteppingAction class
//
//
#include "WLSSteppingAction.hh"

#include "WLSDetectorConstruction.hh"
#include "WLSEventAction.hh"
#include "WLSPhotonDetSD.hh"
#include "WLSSteppingActionMessenger.hh"
#include "WLSUserTrackInformation.hh"

#include "G4ios.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4UImanager.hh"
#include "G4VPhysicalVolume.hh"

// Purpose: Save relevant information into User Track Information

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSSteppingAction::WLSSteppingAction(WLSDetectorConstruction* detector,
                                     WLSEventAction* event)
  : fDetector(detector)
  , fEventAction(event)
{
  fSteppingMessenger = new WLSSteppingActionMessenger(this);
  ResetCounters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSSteppingAction::~WLSSteppingAction() { delete fSteppingMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSSteppingAction::SetBounceLimit(G4int i) { fBounceLimit = i; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::GetNumberOfBounces() { return fCounterBounce; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::GetNumberOfClad1Bounces()
{
  return fCounterClad1Bounce;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::GetNumberOfClad2Bounces()
{
  return fCounterClad2Bounce;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::GetNumberOfWLSBounces() { return fCounterWLSBounce; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::ResetSuccessCounter()
{
  G4int temp  = fCounterEnd;
  fCounterEnd = 0;
  return temp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  G4Track* theTrack = theStep->GetTrack();
  auto trackInformation =
    (WLSUserTrackInformation*) theTrack->GetUserInformation();

  G4StepPoint* thePrePoint  = theStep->GetPreStepPoint();
  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();

  G4VPhysicalVolume* thePrePV  = thePrePoint->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  G4String thePrePVname  = " ";
  G4String thePostPVname = " ";

  if(thePostPV)
  {
    thePrePVname  = thePrePV->GetName();
    thePostPVname = thePostPV->GetName();
  }

  // Recording data for start
  // static const G4ThreeVector ZHat = G4ThreeVector(0.0,0.0,1.0);
  if(theTrack->GetParentID() == 0)
  {
    // This is a primary track
    if(theTrack->GetCurrentStepNumber() == 1)
    {
      //        G4double x  = theTrack->GetVertexPosition().x();
      //        G4double y  = theTrack->GetVertexPosition().y();
      //        G4double z  = theTrack->GetVertexPosition().z();
      //        G4double pz = theTrack->GetVertexMomentumDirection().z();
      //        G4double fInitTheta =
      //                         theTrack->GetVertexMomentumDirection().angle(ZHat);
    }
  }

  // Retrieve the status of the photon
  G4OpBoundaryProcessStatus theStatus = Undefined;

  static G4ThreadLocal G4ProcessManager* OpManager =
    G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if(OpManager)
  {
    G4int nproc = OpManager->GetPostStepProcessVector()->entries();
    G4ProcessVector* fPostStepDoItVector =
      OpManager->GetPostStepProcessVector(typeDoIt);

    for(G4int i = 0; i < nproc; ++i)
    {
      G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
      fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
      if(fOpProcess)
      {
        theStatus = fOpProcess->GetStatus();
        break;
      }
    }
  }

  // Find the skewness of the ray at first change of boundary
  if(fInitGamma == -1 &&
     (theStatus == TotalInternalReflection || theStatus == FresnelReflection ||
      theStatus == FresnelRefraction) &&
     trackInformation->IsStatus(InsideOfFiber))
  {
    G4double px = theTrack->GetVertexMomentumDirection().x();
    G4double py = theTrack->GetVertexMomentumDirection().y();
    G4double x  = theTrack->GetPosition().x();
    G4double y  = theTrack->GetPosition().y();

    fInitGamma = x * px + y * py;

    fInitGamma =
      fInitGamma / std::sqrt(px * px + py * py) / std::sqrt(x * x + y * y);

    fInitGamma = std::acos(fInitGamma * rad);

    if(fInitGamma / deg > 90.0)
    {
      fInitGamma = 180 * deg - fInitGamma;
    }
  }
  // Record Photons that missed the photon detector but escaped from readout
  if(!thePostPV && trackInformation->IsStatus(EscapedFromReadOut))
  {
    //G4cout << "SteppingAction: status = EscapedFromReadOut" << G4endl;
    fEventAction->AddEscaped();
    // UpdateHistogramSuccess(thePostPoint,theTrack);
    ResetCounters();

    return;
  }

  // Assumed photons are originated at the fiber OR
  // the fiber is the first material the photon hits
  switch(theStatus)
  {
    // Exiting the fiber
    case FresnelRefraction:
    case SameMaterial:
      fEventAction->AddExiting();

      if(thePostPVname == "WLSFiber" || thePostPVname == "Clad1" ||
         thePostPVname == "Clad2")
      {
        if(trackInformation->IsStatus(OutsideOfFiber))
          trackInformation->AddStatusFlag(InsideOfFiber);

        // Set the Exit flag when the photon refracted out of the fiber
      }
      else if(trackInformation->IsStatus(InsideOfFiber))
      {
        // EscapedFromReadOut if the z position is the same as fiber's end
        if(theTrack->GetPosition().z() == fDetector->GetWLSFiberEnd())
        {
          trackInformation->AddStatusFlag(EscapedFromReadOut);
          fCounterEnd++;
          fEventAction->AddEscapedEnd();
        }
        else  // Escaped from side
        {
          trackInformation->AddStatusFlag(EscapedFromSide);
          trackInformation->SetExitPosition(theTrack->GetPosition());
          //  UpdateHistogramEscape(thePostPoint,theTrack);

          fCounterMid++;
          fEventAction->AddEscapedMid();
          ResetCounters();
        }

        trackInformation->AddStatusFlag(OutsideOfFiber);
        trackInformation->SetExitPosition(theTrack->GetPosition());
      }

      return;

    // Internal Reflections
    case TotalInternalReflection:

      fEventAction->AddTIR();

      // Kill the track if it's number of bounces exceeded the limit
      if(fBounceLimit > 0 && fCounterBounce >= fBounceLimit)
      {
        theTrack->SetTrackStatus(fStopAndKill);
        trackInformation->AddStatusFlag(murderee);
        ResetCounters();
        G4cout << "\n Bounce Limit Exceeded" << G4endl;
        return;
      }
      break;

    case FresnelReflection:

      fCounterBounce++;
      fEventAction->AddBounce();

      if(thePrePVname == "WLSFiber")
      {
        fCounterWLSBounce++;
        fEventAction->AddWLSBounce();
      }
      else if(thePrePVname == "Clad1")
      {
        fCounterClad1Bounce++;
        fEventAction->AddClad1Bounce();
      }
      else if(thePrePVname == "Clad2")
      {
        fCounterClad2Bounce++;
        fEventAction->AddClad1Bounce();
      }

      // Determine if the photon has reflected off the read-out end
      if(theTrack->GetPosition().z() == fDetector->GetWLSFiberEnd())
      {
        if(!trackInformation->IsStatus(ReflectedAtReadOut) &&
           trackInformation->IsStatus(InsideOfFiber))
        {
          trackInformation->AddStatusFlag(ReflectedAtReadOut);

          if(fDetector->IsPerfectFiber() &&
             theStatus == TotalInternalReflection)
          {
            theTrack->SetTrackStatus(fStopAndKill);
            trackInformation->AddStatusFlag(murderee);
            // UpdateHistogramReflect(thePostPoint,theTrack);
            ResetCounters();
            return;
          }
        }
      }
      return;

    // Reflection off the mirror
    case LambertianReflection:
    case LobeReflection:
    case SpikeReflection:

      fEventAction->AddReflected();
      // Check if it hits the mirror
      if(thePostPVname == "Mirror")
      {
        trackInformation->AddStatusFlag(ReflectedAtMirror);
        fEventAction->AddMirror();
      }
      return;

    // Detected by a detector
    case Detection:
      // Detected automatically with G4OpBoundaryProcess->InvokeSD set true

      // Stop Tracking when it hits the detector's surface
      ResetCounters();
      theTrack->SetTrackStatus(fStopAndKill);
      return;

    default:
      break;
  }

  // Check for absorbed photons
  if(theTrack->GetTrackStatus() != fAlive &&
     trackInformation->IsStatus(InsideOfFiber))
  {
    // UpdateHistogramAbsorb(thePostPoint,theTrack);
    ResetCounters();
    return;
  }
}
