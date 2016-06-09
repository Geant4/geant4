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
// $Id: G4QDiscProcessMixer.hh,v 1.3 2008/07/09 19:45:09 dennis Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//      ---------------- G4QDiscProcessMixer header ----------------
//                 by Mikhail Kossov, Aug 2007.
//  Header of G4QDiscProcessMixer class (hadron+A) of the CHIPS Simulation Branch in GEANT4
// -------------------------------------------------------------------------------
// This is a unique CHIPS class for the Hadron-Nuclear Diffractive Interaction Prosesses
// -------------------------------------------------------------------------------
// At present (Aug-07) it is not tested.
// The normalization is based on the temporary G4QIonIonCrossSection class
// -------------------------------------------------------------------------------

#ifndef G4QDiscProcessMixer_hh
#define G4QDiscProcessMixer_hh

// GEANT4 Headers
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
//#include "G4ParticleTypes.hh"
//#include "G4VParticleChange.hh"
//#include "G4ParticleDefinition.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"
#include "G4QDiscreteProcessVector.hh"
//#include "G4NucleiPropertiesTable.hh"
//#include "G4ThreeVector.hh"
//#include "G4LorentzVector.hh"

#include <vector>

class G4QDiscProcessMixer : public G4VDiscreteProcess
{
public:

  // Constructor
  G4QDiscProcessMixer(const G4String& processName = "Mixed Discrete Process",
                      const G4ParticleDefinition* proj = G4Gamma::Gamma(),
                      G4ProcessType pType = fHadronic );

  // Destructor
  ~G4QDiscProcessMixer();

  G4bool IsApplicable(const G4ParticleDefinition& particle);

  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
			                                             G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep); 

  void AddDiscreteProcess(G4VDiscreteProcess* DP, G4double ME);

  //G4LorentzVector GetEnegryMomentumConservation();

  //G4int GetNumberOfNeutronsInTarget();

private:

  // Hide assignment operator as private 
  G4QDiscProcessMixer& operator=(const G4QDiscProcessMixer &right);

  // Copy constructor
  G4QDiscProcessMixer(const G4QDiscProcessMixer& DPM);

		// BODY
  const G4ParticleDefinition* DPParticle;             // Particle for DiscreteProc mixture
  G4QDiscreteProcessVector theDPVector;               // Vector of Discrete Processes
  //G4LorentzVector EnMomConservation;                // Residual of Energy/Momentum Cons.
  //G4int nOfNeutrons;                                // #of neutrons in the target nucleus
};
#endif
