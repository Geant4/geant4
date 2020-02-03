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
/// \file GB06BOptnSplitAndKillByImportance.hh
/// \brief Definition of the GB06BOptnSplitAndKillByImportance class

#ifndef GB06BOptnSplitAndKillByImportance_hh
#define GB06BOptnSplitAndKillByImportance_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChangeForNothing.hh"
#include "G4TouchableHistoryHandle.hh"
class G4BiasingProcessSharedData;
#include <map>


class GB06BOptnSplitAndKillByImportance : public G4VBiasingOperation {
public:
  // -- Constructor :
  GB06BOptnSplitAndKillByImportance(G4String name);
  // -- destructor:
  virtual ~GB06BOptnSplitAndKillByImportance();

public:
  // ----------------------------------------------
  // -- Methods from G4VBiasingOperation interface:
  // ----------------------------------------------
  // -- Unused:
  virtual const G4VBiasingInteractionLaw*
  ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*,
                                         G4ForceCondition&                 ) final
  {return 0;}
  virtual G4VParticleChange*                            
  ApplyFinalStateBiasing               ( const G4BiasingProcessInterface*,
                                         const G4Track*,
                                         const G4Step*,
                                         G4bool&                           ) final
  {return 0;}

  // -- Used methods ("non-physics biasing methods"):
  // ------------------------------------------------
  // -- Method to return the distance or the condition under which
  // -- requesting the biasing.
  // -- Here we use the condition "forced" and the distance returned
  // -- is "dummy" (DBL_MAX).
  virtual G4double
  DistanceToApplyOperation              ( const G4Track*,
                                          G4double,
                                          G4ForceCondition* condition      ) final;
  // -- Method the generate the final state, which is:
  // --  - made of the primary with half of its original weight, and a clone of it in
  // --    case of splitting
  // --  - the primary with increased weight or the primary killed, in case of killing
  virtual G4VParticleChange* 
  GenerateBiasingFinalState             ( const G4Track*,
                                          const G4Step*                    ) final;

  // -- Specific to this example:
  // ----------------------------
  // -- The parallel world index refers to the index in the list of parallel worlds
  // -- handled by the G4ParallelGeometriesLimiterProcess intance:
  void  SetParallelWorldIndex( G4int parallelWorldIndex )
  { fParallelWorldIndex = parallelWorldIndex; }
  G4int GetParallelWorldIndex() const
  { return fParallelWorldIndex; }
  // -- Method to set the biasing shared data, they contain many
  // -- of the information we need related to parallel geometry:
  void SetBiasingSharedData( const G4BiasingProcessSharedData* sharedData )
  { fBiasingSharedData = sharedData; }
  // -- Set the importance map to be used:
  void SetImportanceMap( std::map< G4int, G4int >* importanceMap )
  { fImportanceMap = importanceMap; }
  
  
private:
  G4int                                   fParallelWorldIndex;
  const G4BiasingProcessSharedData*        fBiasingSharedData;
  G4TouchableHistoryHandle           fPreStepTouchableHistory;
  G4TouchableHistoryHandle          fPostStepTouchableHistory;
  G4ParticleChange                            fParticleChange;
  G4ParticleChangeForNothing             fDummyParticleChange;
  std::map< G4int, G4int >*                    fImportanceMap;
  

};

#endif
