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
// G4AdjointStackingAction
//
// Class description:
//
// Stacking action used in the adjoint simulation. It is responsible to kill
// a primary forward particle before it is traked in the forward phase,
// if the last adjoint particle did not reach the adjoint surface.
// Was needed for the new design where the G4AdjointSimManager is no more an
// extension of G4RunManager. If the primary particles are not killed before
// being tracked in the sensitive geometry, the User Stacking action can be
// used during the forward phase if specified by the method
// G4AdjointSimManager::UseUserStackingAction(G4bool). 

// Author: L. Desorgher, SpaceIT GmbH - April 2008
// Contract: ESA contract 21435/08/NL/AT
// Customer: ESA/ESTEC
// --------------------------------------------------------------------
#ifndef G4AdjointStackingAction_hh
#define G4AdjointStackingAction_hh 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;
class G4ParticleDefinition;
class G4AdjointTrackingAction;

class G4AdjointStackingAction : public G4UserStackingAction
{
  public:

    explicit G4AdjointStackingAction(G4AdjointTrackingAction* anAction);
    ~G4AdjointStackingAction() override = default;

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack) override;
    void NewStage() override;
    void PrepareNewEvent() override;
    inline void SetUserFwdStackingAction(G4UserStackingAction* anAction)
      { theFwdStackingAction = anAction; }
    inline void SetUserAdjointStackingAction(G4UserStackingAction* anAction)
      { theUserAdjointStackingAction = anAction; }
    inline void SetKillTracks(G4bool aBool)
      { kill_tracks =aBool; }
    inline void SetAdjointMode(G4bool aBool)
      { adjoint_mode=aBool; }

  private:

    G4UserStackingAction* theFwdStackingAction = nullptr;
    G4UserStackingAction* theUserAdjointStackingAction = nullptr;
    G4bool reclassification_stage = false;
    G4bool first_reclassification_stage = false;
    G4bool kill_tracks = false;
    G4bool adjoint_mode = false;
    G4AdjointTrackingAction* theAdjointTrackingAction = nullptr;
};

#endif

