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
// $Id: G4ITTrackingManager.cc 103042 2017-03-10 11:50:07Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITTrackingManager.hh"
#include "G4Track.hh"
#include "G4ProcessManager.hh"
#include "G4IT.hh"
#include "G4TrackingInformation.hh"
#include "G4ITTrackingInteractivity.hh"
#include "G4VITSteppingVerbose.hh"
#include "G4ITTrackHolder.hh"

G4ITTrackingManager::G4ITTrackingManager()
{
  fpTrackingInteractivity = 0;
}
//___________________________________________________
//void G4ITTrackingManager::Initialize()
//{
//    if(fpTrackingInteractivity) fpTrackingInteractivity->Initialize();
//}
//___________________________________________________
G4ITTrackingManager::~G4ITTrackingManager()
{
  if (fpTrackingInteractivity) delete fpTrackingInteractivity;
}
//___________________________________________________
void G4ITTrackingManager::StartTracking(G4Track* track)
{
  if (fpTrackingInteractivity)
  {
    fpTrackingInteractivity->StartTracking(track);
#ifdef G4VERBOSE
    fpTrackingInteractivity->GetSteppingVerbose()->TrackingStarted(track);
#endif
  }

  // Inform beginning of tracking to physics processes
  track->GetDefinition()->GetProcessManager()->StartTracking(track);
}
//___________________________________________________
void G4ITTrackingManager::AppendStep(G4Track* track, G4Step* step)
{
  if (fpTrackingInteractivity) fpTrackingInteractivity->AppendStep(track, step);
}

//___________________________________________________
void G4ITTrackingManager::SetInteractivity(G4ITTrackingInteractivity* iteractivity)
{
  if (fpTrackingInteractivity && fpTrackingInteractivity != iteractivity)
  {
    delete fpTrackingInteractivity;
  }
  fpTrackingInteractivity = iteractivity;
}

//___________________________________________________
void G4ITTrackingManager::EndTracking(G4Track* track)
{
  if (fpTrackingInteractivity)
  {
    fpTrackingInteractivity->EndTracking(track);
#ifdef G4VERBOSE
    fpTrackingInteractivity->GetSteppingVerbose()->TrackingEnded(track);
#endif
  }

  G4ITTrackHolder::Instance()->PushToKill(track);
}

void G4ITTrackingManager::EndTrackingWOKill(G4Track* track)
{
  if (fpTrackingInteractivity)
  {
    fpTrackingInteractivity->EndTracking(track);
#ifdef G4VERBOSE
    fpTrackingInteractivity->GetSteppingVerbose()->TrackingEnded(track);
#endif
  }
}
