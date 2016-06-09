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
// $Id: G4ITTrackingManager.hh 65022 2012-11-12 16:43:12Z gcosmo $
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

class G4ITTrackingManager;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4Step;
class G4Track;
class G4ITTrackingInteractivity;

class G4ITTrackingManager
{
protected:
    G4ITTrackingInteractivity* fpTrackingInteractivity;

public:
    G4ITTrackingManager();
    virtual ~G4ITTrackingManager();

    //void Initialize();

    virtual void StartTracking(G4Track*);
    virtual void AppendStep(G4Track* track, G4Step* step);
    virtual void EndTracking(G4Track*);

    void SetInteractivity(G4ITTrackingInteractivity*);
};

#endif // G4ITTRACKINGMANAGER_HH
