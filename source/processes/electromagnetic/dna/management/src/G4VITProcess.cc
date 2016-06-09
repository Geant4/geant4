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
// $Id: G4VITProcess.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
#include "G4VITProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4IT.hh"

size_t G4VITProcess::fNbProcess = 0;

G4VITProcess::G4VITProcess(const G4String& name, G4ProcessType type) :
    G4VProcess( name, type ),
    fpState (0),
    fProcessID(fNbProcess)
{
        fNbProcess++;
        SetInstantiateProcessState(true);
        currentInteractionLength            = 0;
        theInteractionTimeLeft              = 0;
        theNumberOfInteractionLengthLeft    = 0;
        fProposesTimeStep = false;
}

G4VITProcess::G4ProcessState::G4ProcessState()
{
    theNumberOfInteractionLengthLeft = -1.0 ;
    theInteractionTimeLeft           = -1.0 ;
    currentInteractionLength         = -1.0 ;
}

G4VITProcess::G4ProcessState::~G4ProcessState()
{;}

G4VITProcess::~G4VITProcess()
{
    //dtor
    // As the owner, G4IT should delete fProcessState
}

G4VITProcess::G4VITProcess(const G4VITProcess& other)  : G4VProcess(other), fProcessID(other.fProcessID)
{
	//copy ctor
        fpState                             = 0 ;
        currentInteractionLength            = 0;
        theInteractionTimeLeft              = 0;
        theNumberOfInteractionLengthLeft    = 0;
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
        fpState = new G4ProcessState();
    }

    theNumberOfInteractionLengthLeft    = &(fpState->theNumberOfInteractionLengthLeft );
    theInteractionTimeLeft              = &(fpState->theInteractionTimeLeft           );
    currentInteractionLength            = &(fpState->currentInteractionLength         );
    trackingInfo->RecordProcessState(fpState,fProcessID);
    fpState = 0;
}

void G4VITProcess::SubtractNumberOfInteractionLengthLeft(
                                  G4double previousStepSize )
{
  if (fpState->currentInteractionLength>0.0) {
    fpState->theNumberOfInteractionLengthLeft -= previousStepSize/(fpState->currentInteractionLength);
    if(fpState->theNumberOfInteractionLengthLeft<0.) {
       fpState->theNumberOfInteractionLengthLeft=perMillion;
    }

  } else {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VProcess::SubtractNumberOfInteractionLengthLeft()";
      G4cerr << " [" << theProcessName << "]" <<G4endl;
      G4cerr << " currentInteractionLength = " << *currentInteractionLength/cm << " [cm]";
      G4cerr << " previousStepSize = " << previousStepSize/cm << " [cm]";
      G4cerr << G4endl;
    }
#endif
    G4Exception("G4VProcess::SubtractNumberOfInteractionLengthLeft()",
		"Negative currentInteractionLength",EventMustBeAborted,theProcessName);
  }
}
