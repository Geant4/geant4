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
// $Id: G4QDiffraction.hh,v 1.1 2009-11-17 10:36:54 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QDiffraction header ----------------
//                 by Mikhail Kossov, Aug 2007.
//  Header of G4QDiffraction class (hadron+A) of the CHIPS Simulation Branch in GEANT4
// -------------------------------------------------------------------------------
// This is a unique CHIPS class for the Hadron-Nuclear Diffractive Interaction Prosesses
// -------------------------------------------------------------------------------
// At present (Aug-07) it is based on the G4QDiffractionRatio class and is not tested.
// The normalization is based on the temporary G4QProtonNuclearCrossSection class
// -------------------------------------------------------------------------------
// Short description: This is a process, which describes the diffraction
// excitation of the projectile and the nucleus. On nuclei in addition there
// can be a coherent diffraction process for the projectile, but it is
// comparably small. The most important part of the diffraction is the
// progectile diffraction excitation, as in this interaction proton can lose
// only a small part of its energy and make the shower longer. This is because
// only 1-2 (n) pions are produce in the diffraction escitation, and the mean
// kept energy of the nucleon is (1-n/7)=80%. For kaons the kept energy is much
// smaller (1-n/3.5)=60%, and for pions it is less important (about 40%).
// ----------------------------------------------------------------------------

#ifndef G4QDiffraction_hh
#define G4QDiffraction_hh

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
#include "G4QDiffractionRatio.hh"
#include "G4QProtonNuclearCrossSection.hh"
#include "G4QIsotope.hh"
#include "G4QCHIPSWorld.hh"
#include "G4QHadronVector.hh"
#include <vector>

class G4QDiffraction : public G4VDiscreteProcess
{
public:

  // Constructor
  G4QDiffraction(const G4String& processName ="CHIPS_DiffractiveInteraction");

  // Destructor
  ~G4QDiffraction();

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


  G4LorentzVector GetEnegryMomentumConservation();

  G4int GetNumberOfNeutronsInTarget();

private:

  // Hide assignment operator as private 
  G4QDiffraction& operator=(const G4QDiffraction &right);

  // Copy constructor
  G4QDiffraction(const G4QDiffraction&);

  // Calculate Cross-Section of the Diffraction Reaction (p is in GeV @@ units)
  G4double CalculateXS(G4double p, G4int Z, G4int N, G4int pPDG);

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
