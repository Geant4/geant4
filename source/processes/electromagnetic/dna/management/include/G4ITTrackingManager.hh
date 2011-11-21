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

#ifndef G4ITTRACKINGMANAGER_HH
#define G4ITTRACKINGMANAGER_HH

#include "G4TrackList.hh"
/***
#include "G4VTrajectory.hh"
***/

class G4ITTrackingManager;
class G4UserTrackingAction;
class G4Step;
class G4Track;

class G4Trajectory_Lock
{
    friend class G4ITTrackingManager;

    G4Trajectory_Lock()
    /*** : fpTrajectory(0)***/
    {;}
    ~G4Trajectory_Lock()
    {;}

    /***G4VTrajectory* fpTrajectory;***/
};

class G4ITTrackingManager
{
    int fStoreTrajectory;
    int fVerboseLevel;
/***    G4UserTrackingAction* fpUserTrackingAction; ***/
    bool fSetNewUserTrackingAction;

public:
    G4ITTrackingManager();
    ~G4ITTrackingManager();

    void Initialize();

    void StartTracking(G4Track*);
    void TrackBanner(G4Track* track, const G4String& message = "");
    void AppendTrajectory(G4Track* track, G4Step* step);
    void EndTracking(G4Track*);
    inline void SetVerbose(int);

    // if not set, will use tracking action of standard G4TrackingManager
/***
    void SetUserTrackingAction(G4UserTrackingAction*);
    inline G4UserTrackingAction* GetUserTrackingAction();
***/
};

/***
  inline G4UserTrackingAction* G4ITTrackingManager::GetUserTrackingAction()
  {
      return fpUserTrackingAction;
  }
***/

inline void G4ITTrackingManager::SetVerbose(int flag)
{
    fVerboseLevel = flag;
}

#endif // G4ITTRACKINGMANAGER_HH
