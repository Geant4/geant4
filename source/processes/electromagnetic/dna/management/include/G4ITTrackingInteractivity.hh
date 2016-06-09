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
// $Id: G4ITTrackingInteractivity.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
//
// History:
// -----------
//
// -------------------------------------------------------------------

#ifndef G4ITTRACKINGINTERACTIVITY_HH
#define G4ITTRACKINGINTERACTIVITY_HH

#include "G4String.hh"

class G4Track;
class G4Step;
class G4UserTrackingAction;
class G4UserSteppingAction;

class G4ITTrackingInteractivity
{
protected :
    int fVerboseLevel;

public:
    G4ITTrackingInteractivity(){;}
    virtual ~G4ITTrackingInteractivity(){;}

    virtual void StartTracking(G4Track*){;}
    virtual void AppendStep(G4Track* /*track*/, G4Step* /*step*/){;}
    virtual void EndTracking(G4Track*){;}

    virtual void TrackBanner(G4Track* /*track*/, const G4String& message = "");
    inline void SetVerbose(int);
};

inline void G4ITTrackingInteractivity::SetVerbose(int flag)
{
    fVerboseLevel = flag;
}

#endif // G4ITTRACKINGINTERACTIVITY_HH
