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
// $Id$
//
//      ---------------- G4QIonIonElastic header ----------------
//                 by Mikhail Kossov, December 2006.
//  Header of G4QIonIonElastic class (A+A) of the CHIPS Simulation Branch in GEANT4
// -------------------------------------------------------------------------------
// This is a unique CHIPS class for the Ion-Ion Elastic Scattering Prosesses
// -------------------------------------------------------------------------------
// @@ This class is on the testing level @@
// -------------------------------------------------------------------------------
// Short description: a simple process for the Ion-Ion elastic scattering.
// For heavy by heavy ions it can reach 50% of the total cross-section.
// -----------------------------------------------------------------------

#ifndef G4QIonIonElastic_hh
#define G4QIonIonElastic_hh

// GEANT4 Headers
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleTypes.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

// CHIPS Headers
#include "G4QIonIonCrossSection.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"
#include "G4QIsotope.hh"
#include "G4QCHIPSWorld.hh"
#include "G4QHadron.hh"
#include <vector>

class G4QIonIonElastic : public G4VDiscreteProcess
{
public:

  // Constructor
  G4QIonIonElastic(const G4String& processName ="CHIPS_IonIonElasticScattering");

  // Destructor
  ~G4QIonIonElastic();

  G4bool IsApplicable(const G4ParticleDefinition& particle);

  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                           G4ForceCondition* condition);
  // It returns the MeanFreePath of the process for the current track :
  // (energy, material)
  // The previousStepSize and G4ForceCondition* are not used.
  // This function overloads a virtual function of the base class.        
  // It is invoked by the ProcessManager of the Particle.
 

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep); 
  // It computes the final state of the process (at end of step),
  // returned as a ParticleChange object.       
  // This function overloads a virtual function of the base class.
  // It is invoked by the ProcessManager of the Particle.


  G4LorentzVector GetEnegryMomentumConservation(){return EnMomConservation;}

  G4int GetNumberOfNeutronsInTarget(){return nOfNeutrons;}

private:

  // Hide assignment operator as private 
  G4QIonIonElastic& operator=(const G4QIonIonElastic &right);

  // Copy constructor
  G4QIonIonElastic(const G4QIonIonElastic&);

  // BODY
  // Static Parameters --------------------------------------------------------------------
  static G4int    nPartCWorld; // The#of particles for hadronization (limit of A of fragm.)
  //--------------------------------- End of static parameters ---------------------------
  // Working parameters
  G4VQCrossSection* theCS;
  G4LorentzVector EnMomConservation;                  // Residual of Energy/Momentum Cons.
  G4int nOfNeutrons;                                  // #of neutrons in the target nucleus

  // Modifires for the reaction
  G4double Time;                                      // Time shift of the capture reaction
  G4double EnergyDeposition;                          // Energy deposited in the reaction
  static std::vector <G4int> ElementZ;                // Z of the element(i) in theLastCalc
  static std::vector <G4double> ElProbInMat;          // SumProbabilityElements in Material
  static std::vector <std::vector<G4int>*> ElIsoN;    // N of isotope(j) of Element(i)
  static std::vector <std::vector<G4double>*> IsoProbInEl;// SumProbabIsotopes in Element i
};
#endif

