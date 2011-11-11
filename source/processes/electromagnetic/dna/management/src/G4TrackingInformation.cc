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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4TrackingInformation.hh"
#include "G4VITProcess.hh"

static const size_t SizeOfSelectedDoItVector=100;

G4TrackingInformation::G4TrackingInformation() :
    fStepLeader                     (false),
    fStepStatus                     (fUndefined),
    fPhysicalStep                   (-1.),
    fTouchableHandle                (0),
    fSelectedAtRestDoItVector       ((size_t)G4VITProcess::GetMaxProcessIndex(),0),
    fSelectedAlongStepDoItVector    ((size_t)G4VITProcess::GetMaxProcessIndex(),0),
    fSelectedPostStepDoItVector     ((size_t)G4VITProcess::GetMaxProcessIndex(),0),
    fProcessState                   ((size_t)G4VITProcess::GetMaxProcessIndex(),0)
{
    //ctor
    fpTrajectory_Lock = 0;
    fRecordedTrackGlobalTime = -1;
    fRecordedTrackLocalTime = -1;
    fSafety = -1;
    fPreviousStepSize = -1;

}

G4TrackingInformation::~G4TrackingInformation()
{
    //dtor
    for(int i = 0 ; i < (int) fProcessState.size() - 1 ; i++)
    {
        if(fProcessState[i])
        {
            delete fProcessState[i];
            fProcessState[i] = 0;
        }
    }
    fProcessState.clear();
}

// should not be used
G4TrackingInformation::G4TrackingInformation(const G4TrackingInformation& /*other*/) :
    fStepLeader                     (false),
    fStepStatus                     (fUndefined),
    fPhysicalStep                   (-1.),
    fTouchableHandle                (0)
{
    //copy ctor
    fpTrajectory_Lock = 0;
    fRecordedTrackGlobalTime = -1;
    fRecordedTrackLocalTime = -1;
    fSafety = -1;
    fPreviousStepSize = -1;
}

// should not be used
G4TrackingInformation& G4TrackingInformation::operator=(const G4TrackingInformation& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

G4ProcessState_Lock* G4TrackingInformation::GetProcessState(G4int index)
{
    if(index> G4VITProcess::GetMaxProcessIndex() || index<0)
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "G4TrackingInformation::GetProcInfo : Wrong process subType : " ;
        exceptionDescription << index ;
        G4Exception("G4TrackingInformation::GetProcessState","G4TrackingInformation003",
                    FatalErrorInArgument,exceptionDescription);
    }

    return fProcessState[index];
}

void G4TrackingInformation::RecordCurrentPositionNTime(G4Track* track)
{
    if(track)
    {
        fRecordedTrackPosition = track->GetPosition();
        fRecordedTrackLocalTime = track->GetLocalTime();
        fRecordedTrackGlobalTime = track->GetGlobalTime();
    }
}

