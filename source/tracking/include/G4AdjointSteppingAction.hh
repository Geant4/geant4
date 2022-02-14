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
// G4AdjointSteppingAction
//
// Class description:
//
// Stepping action used in the adjoint simulation. 
// It is responsible to stop the adjoint tracking phase when:
// a) The adjoint track reaches the external surface.  
// b) The being tracked adjoint dynamic particle get an energy higher than
//    the maximum energy of the external source.
// c) The adjoint track enters the volume delimited by the adjoint source.
// In the case a) the info (energy,weight,...) of the adjoint dynamic particle
// associated to the track when crossing the external source is registered and
// in the next event a forward primary is generated.
// In the other cases b) and c) the next generated fwd particle is killed
// before being tracked and the next tracking of an adjoint particle is
// started directly.

// Author: L. Desorgher, SpaceIT GmbH
// Contract: ESA contract 21435/08/NL/AT
// Customer: ESA/ESTEC
// History: 
// - 15/01/2007 L.Desorgher, created.
// - 04/11/2009 L.Desorgher, added possibility to use user stepping action.
// - 20/11/2009 L.Desorgher, corrected stop of adjoint particles tracking
//                           when reentering the adjoint source.
//---------------------------------------------------------------------
#ifndef G4AdjointSteppingAction_hh
#define G4AdjointSteppingAction_hh 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4AdjointCrossSurfChecker;
class G4ParticleDefinition;

class G4AdjointSteppingAction : public G4UserSteppingAction
{
  public:

    G4AdjointSteppingAction();
   ~G4AdjointSteppingAction();

    void UserSteppingAction(const G4Step*);
    
    inline void SetExtSourceEMax(G4double Emax)
           { ext_sourceEMax = Emax; } 
    inline void SetStartEvent(G4bool aBool)
           { start_event = aBool; }
    inline G4bool GetDidAdjParticleReachTheExtSource()
           { return did_adj_part_reach_ext_source; }
    inline G4ThreeVector GetLastMomentum()
           { return last_momentum; }
    inline G4ThreeVector GetLastPosition()
           { return last_pos; }
    inline G4double GetLastEkin()
           { return last_ekin; }
    inline G4double GetLastWeight()
           { return last_weight; }
    inline void SetPrimWeight(G4double weight)
           { prim_weight = weight; } 
    inline G4ParticleDefinition* GetLastPartDef()
           { return last_part_def; }
    inline void SetUserAdjointSteppingAction(G4UserSteppingAction* anAction)
           { theUserAdjointSteppingAction = anAction; }
    inline void SetUserForwardSteppingAction(G4UserSteppingAction* anAction)
           { theUserFwdSteppingAction = anAction; }
    inline void SetAdjointTrackingMode(G4bool aBool)
           { is_adjoint_tracking_mode = aBool; }
    inline void ResetDidOneAdjPartReachExtSourceDuringEvent()
           { did_one_adj_part_reach_ext_source_during_event = false; }
    inline void SetAdjointGeantinoTrackingMode(G4bool aBool)
           { is_adjoint_geantino_tracking_mode = aBool; }

  private:

    G4double ext_sourceEMax = 0.0;
    G4AdjointCrossSurfChecker* theG4AdjointCrossSurfChecker = nullptr;
    
    G4ThreeVector last_momentum, last_pos; 
    G4double last_ekin = 0.0;
    G4double last_weight = 0.0;
    G4double prim_weight = 1.0;
    G4ParticleDefinition* last_part_def = nullptr;
    G4UserSteppingAction* theUserAdjointSteppingAction = nullptr;
    G4UserSteppingAction* theUserFwdSteppingAction = nullptr;

    G4bool start_event = false;
    G4bool did_adj_part_reach_ext_source = false;
    G4bool did_one_adj_part_reach_ext_source_during_event = false;
    G4bool is_adjoint_tracking_mode = false;
    G4bool is_adjoint_geantino_tracking_mode = false;
};

#endif
