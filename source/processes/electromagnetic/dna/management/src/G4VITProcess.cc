#include "G4VITProcess.hh"
#include "G4IT.hh"

G4int G4VITProcess::fNbProcess = 0;

G4VITProcess::G4VITProcess(const G4String& name, G4ProcessType type) :
    G4VProcess( name, type ),
    fState (0),
    fProcessID(fNbProcess)
{
        fNbProcess++;
        SetInstantiateProcessState(true);
        currentInteractionLength            = 0;
        theInteractionTimeLeft              = 0;
        theNumberOfInteractionLengthLeft    = 0;
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
        fState                              = 0 ;
        currentInteractionLength            = 0;
        theInteractionTimeLeft              = 0;
        theNumberOfInteractionLengthLeft    = 0;
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
        fState = new G4ProcessState();
    }

    theNumberOfInteractionLengthLeft    = &(fState->theNumberOfInteractionLengthLeft );
    theInteractionTimeLeft              = &(fState->theInteractionTimeLeft           );
    currentInteractionLength            = &(fState->currentInteractionLength         );
    trackingInfo->RecordProcessState(fState,fProcessID);
    fState = 0;
}

void G4VITProcess::SubtractNumberOfInteractionLengthLeft(
                                  G4double previousStepSize )
{
  if (fState->currentInteractionLength>0.0) {
    fState->theNumberOfInteractionLengthLeft -= previousStepSize/(fState->currentInteractionLength);
    if(fState->theNumberOfInteractionLengthLeft<0.) {
       fState->theNumberOfInteractionLengthLeft=perMillion;
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
