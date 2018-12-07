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
//
/// \file GB05BOptnSplitAndKillByCrossSection.hh
/// \brief Definition of the GB05BOptnSplitAndKillByCrossSection class

#ifndef GB05BOptnSplitAndKillByCrossSection_hh
#define GB05BOptnSplitAndKillByCrossSection_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ParticleChange.hh"


class GB05BOptnSplitAndKillByCrossSection : public G4VBiasingOperation {
public:
  // -- Constructor :
  GB05BOptnSplitAndKillByCrossSection(G4String name);
  // -- destructor:
  virtual ~GB05BOptnSplitAndKillByCrossSection();

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
  // -- Here this distance will be sampled according the exponential
  // -- interaction law, using the interaction length passed to the
  // -- method SetInteractionLength(G4double)  below.
  virtual G4double
  DistanceToApplyOperation              ( const G4Track*,
                                          G4double,
                                          G4ForceCondition* condition      ) final;
  // -- Method the generate the final state, which is made of the primary
  // -- with half of its original weight, and a clone of it.
  virtual G4VParticleChange* 
  GenerateBiasingFinalState             ( const G4Track*,
                                          const G4Step*                    ) final;

  // -- Specific to this example:
  // ----------------------------
  void SetInteractionLength(G4double interactionLength )
  {
    fInteractionLength = interactionLength;
  }

  
private:
  G4ParticleChange    fParticleChange;
  G4double         fInteractionLength;

};

#endif
