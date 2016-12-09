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
// $Id: $
//
//---------------------------------------------------------------
//
// G4BOptnForceFreeFlight
//
// Class Description:
//    A G4VBiasingOperation physics-based biasing operation.
//    If forces the physics process to not act on the track.
//    In this implementation (meant for the ForceCollision
//    operator) the free flight is done under zero weight for
//    the track, and the action is meant to accumulate the weight
//    change for making this uninteracting flight,
//    cumulatedWeightChange.
//    When the track reaches the current volume boundary, its
//    weight is restored with value :
//            initialWeight * cumulatedWeightChange
//
//---------------------------------------------------------------
//   Initial version                         Nov. 2013 M. Verderi
#ifndef G4BOptnForceFreeFlight_hh
#define G4BOptnForceFreeFlight_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ForceCondition.hh"
#include "G4ParticleChange.hh" // -- §§ should add a dedicated "weight change only" particle change
class G4ILawForceFreeFlight;


class G4BOptnForceFreeFlight : public G4VBiasingOperation {
public:
  // -- Constructor :
  G4BOptnForceFreeFlight(G4String name);
  // -- destructor:
  virtual ~G4BOptnForceFreeFlight();
  
public:
  // -- Methods from G4VBiasingOperation interface:
  // -------------------------------------------
  // -- Used:
  virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*, G4ForceCondition& );
  virtual void                                                      AlongMoveBy( const G4BiasingProcessInterface*, const G4Step*, G4double );
  virtual G4VParticleChange*                             ApplyFinalStateBiasing( const G4BiasingProcessInterface*, const G4Track*, const G4Step*, G4bool&);

  // -- Unused:
  virtual G4double                               DistanceToApplyOperation( const G4Track*,
									   G4double,
									   G4ForceCondition*)  {return DBL_MAX;}
  virtual G4VParticleChange*                    GenerateBiasingFinalState( const G4Track*,
									   const G4Step*    )  {return 0;}


public:
  // -- Additional methods, specific to this class:
  // ----------------------------------------------
  // -- return concrete type of interaction law:
  G4ILawForceFreeFlight* GetForceFreeFlightLaw() {
    return fForceFreeFlightInteractionLaw;
  }
  // -- initialization for weight:
  void ResetInitialTrackWeight(G4double w) {fInitialTrackWeight = w; fCumulatedWeightChange = 1.0;}
  G4bool OperationComplete() const { return fOperationComplete; }
  
private:
  G4ILawForceFreeFlight* fForceFreeFlightInteractionLaw;
  G4double               fCumulatedWeightChange,
                         fInitialTrackWeight;
  G4ParticleChange       fParticleChange;
  G4bool                 fOperationComplete;
};

#endif
