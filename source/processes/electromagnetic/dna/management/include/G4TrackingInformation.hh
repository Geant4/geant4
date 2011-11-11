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
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4TRACKINGINFORMATION_HH
#define G4TRACKINGINFORMATION_HH

#include "globals.hh"
#include <vector>
#include "G4StepStatus.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"

class G4ITStepProcessor;

typedef std::vector<G4int>
        G4SelectedAtRestDoItVector;
typedef std::vector<G4int>
        G4SelectedAlongStepDoItVector;
typedef std::vector<G4int>
        G4SelectedPostStepDoItVector;

class G4Trajectory_Lock;
class G4Track;
struct G4ProcessState_Lock;

/** The class G4TrackingInformation (hold by G4IT)
 *  emcompasses processes informations computed
 *  at the PS/AS/AtRest/InteractionLength stage,
 *  and also, the selection of processes for the
 *  given step.
 */
class G4TrackingInformation
{
public:
    G4TrackingInformation();
    ~G4TrackingInformation();

    //________________________________________________
    /** If the track is the one having the minimum step time,
         *  then it "leads" the step. It will interact will all the
         *  other tracks will be transported.
         */
    inline bool IsLeadingStep(){return fStepLeader;}
    inline void SetLeadingStep(bool value ){fStepLeader = value;}

    //________________________________________________
    /** Every process should store the information
         * computed at the InteractionLegth stage in the track.
         */

    G4ProcessState_Lock* GetProcessState(G4int index);

    inline void RecordProcessState(G4ProcessState_Lock*,
                                   G4int index);

    inline G4Trajectory_Lock* GetTrajectory_Lock()
    {
        return fpTrajectory_Lock ;
    }

    inline void SetTrajectory_Lock(G4Trajectory_Lock* trajLock)
    {
        fpTrajectory_Lock = trajLock;
    }

    void RecordCurrentPositionNTime(G4Track*);
    inline const G4ThreeVector&    GetPreStepPosition() const;
    inline G4double         GetPreStepLocalTime() const;
    inline G4double         GetPreStepGlobalTime() const;

    //-------------
protected:
    //-------------
    friend class G4ITStepProcessor;
    //_______________________________________________________
    G4bool                          fStepLeader ;
    G4StepStatus                    fStepStatus;

    //_______________________________________________________
    G4double    fPhysicalStep ;
    G4double    fSafety;
    G4double    fPreviousStepSize;
    //_______________________________________________________
    G4TouchableHandle fTouchableHandle;

    //_______________________________________________________
    G4Trajectory_Lock* fpTrajectory_Lock;

    //_______________________________________________________
    G4ThreeVector   fRecordedTrackPosition;
    G4double        fRecordedTrackLocalTime;
    G4double        fRecordedTrackGlobalTime;

    //_______________________________________________________
    G4SelectedAtRestDoItVector      fSelectedAtRestDoItVector;
    G4SelectedAlongStepDoItVector   fSelectedAlongStepDoItVector;
    G4SelectedPostStepDoItVector    fSelectedPostStepDoItVector;

    //_______________________________________________________
    /** Holds the information related to processes
      *  Indexed on GetPhysIntVector
      * (cf. G4ITStepProcessor header)
      */
    std::vector<G4ProcessState_Lock*> fProcessState;

    //_______________________________________________________
    /** Copy constructor
         *  \param other Object to copy from
         */
    G4TrackingInformation(const G4TrackingInformation& other);

    /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
    G4TrackingInformation& operator=(const G4TrackingInformation& other);
};

inline void G4TrackingInformation::RecordProcessState(G4ProcessState_Lock* state,
                               G4int index)
{
    fProcessState[index] = state;
}


inline G4double G4TrackingInformation::GetPreStepGlobalTime() const
{
    return fRecordedTrackGlobalTime;
}

inline G4double G4TrackingInformation::GetPreStepLocalTime() const
{
    return fRecordedTrackLocalTime;
}

inline const G4ThreeVector& G4TrackingInformation::GetPreStepPosition() const
{
    return fRecordedTrackPosition;
}


#endif // G4TRACKINGINFORMATION_HH
