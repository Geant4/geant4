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
//---------------------------------------------------------------------
//
// G4BOptnLeadingParticle
//
// Class Description:
//       A G4VBiasingOperation that implements the so-called "Leading
//   particle biasing scheme". It is of interest in the shield problem
//   to estimate the flux leaking from the shield.
//       It works as follows:
//       - it is intented for hadronic inelastic interaction
//       - at each interaction, are kept:
//           - the most energetic particle (the leading particle)
//               - with unmodified weight
//           - randomly one particle of each species
//               - with this particle weight = n * primary_weight where
//                 n is the number of particles of this species
//---------------------------------------------------------------------
// Initial version                                 Nov. 2019 M. Verderi


#ifndef G4BOptnLeadingParticle_hh
#define G4BOptnLeadingParticle_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ParticleChange.hh"

class G4BOptnLeadingParticle : public G4VBiasingOperation {
public:
  // -- Constructor :
  G4BOptnLeadingParticle(G4String name);
  // -- destructor:
  virtual ~G4BOptnLeadingParticle();
  
public:
  // -- Methods from G4VBiasingOperation interface:
  // ----------------------------------------------
  // -- Unused:
  virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*, G4ForceCondition& ) {return nullptr;}
  // -- Used:
  virtual G4VParticleChange*                             ApplyFinalStateBiasing( const G4BiasingProcessInterface*, // -- Method used for this biasing. The related biasing operator
										 const G4Track*,                   // -- returns this biasing operation at the post step do it level
										 const G4Step*,                    // -- when the wrapped process has won the interaction length race.
										 G4bool& );                        // -- The wrapped process final state is then trimmed.
  // -- Unused:
  virtual G4double                                     DistanceToApplyOperation( const G4Track*,
										 G4double,
										 G4ForceCondition*)                                    {return 0;}
  virtual G4VParticleChange*                         GenerateBiasingFinalState( const G4Track*,
										const G4Step*  )                                       {return nullptr;}

public:
  // -- The possibility is given to further apply a Russian roulette on tracks that are accompagnying the leading particle
  // -- after the classical leading particle biasing algorithm has been applied.
  // -- This is of interest when applying the technique to e+ -> gamma gamma for example. Given one gamma is leading,
  // -- the second one is alone in its category, hence selected. With the Russian roulette it is then possible to keep
  // -- this one randomly. This is also of interest for pi0 decays, or for brem. e- -> e- gamma where the e- or gamma
  // -- are alone in their category.
  void     SetFurtherKillingProbability( G4double p )       { fRussianRouletteKillingProbability = p;    }   // -- if p <= 0.0 the killing is ignored.
  G4double GetFurtherKillingProbability()             const { return fRussianRouletteKillingProbability; }
  
private:
  // -- Particle change used to return the trimmed final state:
  G4ParticleChange fParticleChange;
  G4double         fRussianRouletteKillingProbability;

  
};

#endif
