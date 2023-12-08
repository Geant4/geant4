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
#include "G4VITProcess.hh"

#include "G4IT.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackingInformation.hh"

#include <memory>

/*G4ThreadLocal*/size_t *G4VITProcess::fNbProcess = nullptr;

G4VITProcess::G4VITProcess(const G4String& name, G4ProcessType type) :
    G4VProcess(name, type)
//    fpState (0)//,
//fProcessID(fNbProcess)
{
  fpState.reset();
  if(fNbProcess == nullptr) fNbProcess = new size_t(0);
  fProcessID = *fNbProcess;
  (*fNbProcess)++;
  SetInstantiateProcessState(true);
  currentInteractionLength = nullptr;
  theInteractionTimeLeft = nullptr;
  theNumberOfInteractionLengthLeft = nullptr;
  fProposesTimeStep = false;
}

G4VITProcess::G4ProcessState::G4ProcessState()
{
  theNumberOfInteractionLengthLeft = -1.0;
  theInteractionTimeLeft = -1.0;
  currentInteractionLength = -1.0;
}

G4VITProcess::G4ProcessState::~G4ProcessState()
= default;

G4VITProcess::~G4VITProcess()
= default;

G4VITProcess::G4VITProcess(const G4VITProcess& other) :
    G4VProcess(other), fProcessID(other.fProcessID)
{
  //copy ctor
  //fpState                             = 0 ;
  currentInteractionLength            = nullptr;
  theInteractionTimeLeft              = nullptr;
  theNumberOfInteractionLengthLeft    = nullptr;
  fInstantiateProcessState            = other.fInstantiateProcessState;
  fProposesTimeStep                   = other.fProposesTimeStep;
}

G4VITProcess& G4VITProcess::operator=(const G4VITProcess& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

void G4VITProcess::StartTracking(G4Track* track)
{
  G4TrackingInformation* trackingInfo = GetIT(track)->GetTrackingInfo();
  if(InstantiateProcessState())
  {
    //        fpState = new G4ProcessState();
    fpState = std::make_shared<G4ProcessState>();
  }

  theNumberOfInteractionLengthLeft    = &(fpState->theNumberOfInteractionLengthLeft );
  theInteractionTimeLeft              = &(fpState->theInteractionTimeLeft           );
  currentInteractionLength            = &(fpState->currentInteractionLength         );
  trackingInfo->RecordProcessState(fpState,fProcessID);
  //    fpState = 0;
  fpState.reset();
}

