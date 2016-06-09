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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QNGamma header ----------------
//                 by Mikhail Kossov, April 2012.
//  Header of G4QNGamma class (neutrons) of the CHIPS Simulation Branch in GEANT4
// -------------------------------------------------------------------------------
// This is a CHIPS class for the (n,gamma) reactions.
// -------------------------------------------------------------------------------
// ****************************************************************************************
// This Header is a part of the CHIPS physics package (author: M. Kosov)
// ****************************************************************************************
// Short description: This is a CHIPS process class for (n,gamma) reactions.
// -------------------------------------------------------------------------

#ifndef G4QNGamma_hh
#define G4QNGamma_hh

// GEANT4 Headers
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh" 
#include "G4QThd.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleTypes.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4RandomDirection.hh"

// CHIPS Headers
#include "G4VQCrossSection.hh"
#include "G4QIsotope.hh"
#include "G4QNeutronNuclearCrossSection.hh"
#include "G4QNeutronCaptureRatio.hh"
#include "G4QPDGToG4Particle.hh"

class G4QNGamma : public G4VDiscreteProcess
{
public:

  // Constructor
  G4QNGamma(const G4String& processName = "CHIPS_N-Gamma" );

  // Destructor
  ~G4QNGamma();

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

  // Fake void functions
  void SetPhysicsTableBining(G4double, G4double, G4int) {;}
  void BuildPhysicsTable(const G4ParticleDefinition&) {;}
  void PrintInfoDefinition() {;}

  // Internal Energy-Momentum Residual
  G4LorentzVector GetEnegryMomentumConservation();

  // Number of neutrons in the target nucleus (primary)
  G4int GetNumberOfNeutronsInTarget();


private:

  // Hide assignment operator as private 
  G4QNGamma& operator=(const G4QNGamma &right);

  // Copy constructor
  G4QNGamma(const G4QNGamma&);

  // BODY
  // Working parameters
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
