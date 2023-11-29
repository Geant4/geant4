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
/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SteppingAction_H
#define SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction {
  public:
    SteppingAction();
    virtual ~SteppingAction();

    // This is the main method where the properties of the primary particle,
    // at the beginning and when it decays, are collected or computed, and
    // then filled up in the Run object where they are stored (and then
    // printed out at the end of the Run).
    // (For simplicity and brevity, we avoid histograms and compute instead
    //  some statistics ourself, which will be print-out at the end of the run.)
    virtual void UserSteppingAction( const G4Step* ) override;
  
    // This method is called by RunAction::BeginOfRunAction for the
    // initialization of the stepping-action at the beginning of each Run.
    // This is necessary because different runs can have different primary particle
    // types, kinetic energies, starting positions and/or directions, as well as
    // different detector configurations.
    void Initialize();

    // This method is called by RunAction::BeginOfRunAction for providing to the
    // stepping-action the pointer to the run object at the beginning of each Run.
    // This pointer is then used to provide information on the primary decays
    // by calling accessor methods of the Run object.
    void SetRunPointer( Run* inputValue = nullptr ) { fRunPtr = inputValue; }
     
  private:

    G4bool IsPreassignedDecayEnabled() const { return fIsPreassignedDecayEnabled; }
    G4bool IsBoostToLabEnabled() const { return fIsBoostToLabEnabled; }
    G4double ToleranceEPviolations() const { return fToleranceEPviolations; }
    G4double ToleranceDeltaDecayRadius() const { return fToleranceDeltaDecayRadius; }
  
    // Decide whether you want to preassign a decay to the primary particle, if it is unstable.
    // Normally, the preassigned decay is set by a Monte Carlo Event Generator.
    // Here, instead, just for testing, we use the Geant4 decay table, selecting randomly
    // one of the decay channels defined in Geant4 for this primary particle.
    // Here we have two choices: either to set the decay products directly as returned from the
    // decay channel, i.e in the parent particle rest frame; or to boost them in lab frame.
    // Both cases should provide the same result, because internally Geant4 knows whether
    // the reference frame is at test ( E_total = mass of the parent, decaying particle)
    // or the laboratory frame.
    const G4bool fIsPreassignedDecayEnabled = true;  //***LOOKHERE***
    const G4bool fIsBoostToLabEnabled = true;        //***LOOKHERE***

    // Energy-momentum violations are tolerated if smaller than this value
    G4double fToleranceEPviolations;
  
    // Differences between the MC-true decay radius and the real decay radius are tolerated
    // if smaller than this value 
    const G4double fToleranceDeltaDecayRadius = 1.0*CLHEP::micrometer;

    Run* fRunPtr;  // Pointer to the Run object

    // Information regarding the primary particle at the beginning of its tracking 
    G4int    fPrimaryParticleId;
    G4double fPrimaryParticleInitialKineticEnergy;
    G4double fPrimaryParticleInitialTotalEnergy;
    G4double fPrimaryParticleInitialMomentum;
    G4double fPrimaryParticleInitialBeta;     // Lorentz beta
    G4double fPrimaryParticleInitialGamma;    // Lorentz gamma
    G4ThreeVector fPrimaryParticleInitial3Momentum;
    G4ThreeVector fPrimaryParticleInitialPosition;

    // Discrepancies between alternative, but equivalent ways to compute kinematical properties
    // of the primary particle, due to limited numerical accuracy. 
    G4double fMaxEkin_deltaMax;            
    G4double fMaxEtot_deltaMax;            
    G4double fMaxP_deltaMax;               
    G4double fMaxPdir_deltaMax;
    G4double fMaxMass_deltaMax1;           
    G4double fMaxMass_deltaMax2;           
    G4double fMaxMass_deltaMax3;           
    G4double fMeanMass_deltaMax3;          
    G4double fMaxBeta_deltaMax1;           
    G4double fMaxBeta_deltaMax2;           
    G4double fMaxGamma_deltaMax1;          
    G4double fMaxGamma_deltaMax2;          
    G4double fMaxGamma_deltaMax3;          
    G4double fMaxT_proper_deltaMax;        
    G4double fMaxT_lab_deltaMax ;          
    G4double fMaxMc_truth_rPos_deltaMax;   
    G4double fMeanMc_truth_rPos_deltaMax;

    // Properties of the primary particle at the moment of its decay 
    G4double fMeanDeltaR_primaryDecay;   
    G4double fMinDeltaR_primaryDecay;    
    G4double fMaxDeltaR_primaryDecay;    
    G4double fMeanR_primaryDecay;        
    G4double fMinR_primaryDecay;         
    G4double fMaxR_primaryDecay;         
    G4double fMeanX_primaryDecay;        
    G4double fMinX_primaryDecay;         
    G4double fMaxX_primaryDecay;         
    G4double fMeanY_primaryDecay;        
    G4double fMinY_primaryDecay;         
    G4double fMaxY_primaryDecay;         
    G4double fMeanZ_primaryDecay;        
    G4double fMinZ_primaryDecay;         
    G4double fMaxZ_primaryDecay;
    G4double fMeanDeltaAngle_primaryDecay;
    G4double fMinDeltaAngle_primaryDecay;
    G4double fMaxDeltaAngle_primaryDecay;
    G4double fMeanDeltaEkin_primaryDecay;
    G4double fMinDeltaEkin_primaryDecay; 
    G4double fMaxDeltaEkin_primaryDecay; 
    G4double fMeanEkin_primaryDecay;     
    G4double fMinEkin_primaryDecay;      
    G4double fMaxEkin_primaryDecay;      
    G4double fMeanPx_primaryDecay;       
    G4double fMinPx_primaryDecay;        
    G4double fMaxPx_primaryDecay;        
    G4double fMeanPy_primaryDecay;       
    G4double fMinPy_primaryDecay;        
    G4double fMaxPy_primaryDecay;        
    G4double fMeanPz_primaryDecay;       
    G4double fMinPz_primaryDecay;        
    G4double fMaxPz_primaryDecay; 

    // Conceptually wrong ways to estimate the "MC-truth" decay radius and the decay radius
    // (nevertheless useful to estimate the expected potential errors that users might do)
    G4double fMinUnderestimated_mc_truth_rPos_delta;  
    G4double fMaxOverestimated_mc_truth_rPos_delta;   
    G4double fMeanUnderestimated_mc_truth_rPos_delta; 
    G4double fMeanOverestimated_mc_truth_rPos_delta;  
    G4double fMinUnderestimated_rDeltaPos;            
    G4double fMaxOverestimated_rDeltaPos;             
    G4double fMeanUnderestimated_rDeltaPos;           
    G4double fMeanOverestimated_rDeltaPos;            
 
    // Max error due to the use of  float  instead of  double
    G4double fMaxFloat_rDeltaPos_deltaMax;

    // Energy-momentum violation in the decay of the primary particle, computed as difference
    // between the sum of its daughters and the parent (at the moment of its decay)
    G4double fMeanViolationE_primaryDecay; 
    G4double fMinViolationE_primaryDecay;  
    G4double fMaxViolationE_primaryDecay;  
    G4double fMeanViolationPx_primaryDecay;
    G4double fMinViolationPx_primaryDecay; 
    G4double fMaxViolationPx_primaryDecay; 
    G4double fMeanViolationPy_primaryDecay;
    G4double fMinViolationPy_primaryDecay; 
    G4double fMaxViolationPy_primaryDecay; 
    G4double fMeanViolationPz_primaryDecay;
    G4double fMinViolationPz_primaryDecay; 
    G4double fMaxViolationPz_primaryDecay; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
