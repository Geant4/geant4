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
/// \file Run.hh
/// \brief Definition of the Run class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run {
  // This class accumulates relevant quantities (related to the primay decays)
  // collected during the run.
  // ( Note: these information are provided via calls of accessor methods of
  //         this Run class made by SteppingAction::UserSteppingAction. )
  // At the end of a run, the  printInfo  method is called by the run-action
  // to print out some summary information (typically: mean, min and max) about
  // these quantities. In multithreaded (MT) mode, an object of this class is
  // filled up for each working thread, and then merged (automatically by the
  // Geant4 kernel) into another object (of this class) owned by the master class;
  // the  printInfo  method is then called only for the latter run object.
  // Note that, for simplicity and brevity, we avoid histograms and print-out
  // instead some statistics (compute by ourself) at the end of the run.
  public:
  
    Run();
    ~Run();
  
    // This method is called automatically by the Geant4 kernel (not by the user!)
    // at the end of each event. In the case of multithreaded mode, it is called
    // only for the Working thread that handled that event.
    virtual void RecordEvent( const G4Event* anEvent ) override;

    // This method is called automatically by the Geant4 kernel (not by the user!)
    // only in the case of multithreaded mode and only for Working threads.
    virtual void Merge( const G4Run* aRun ) override;
  
    // This method is called by RunAction::EndOfRunAction : in the case
    // of multithreaded mode, only the master thread calls it.
    void PrintInfo() const;

    void SetPrimaryParticleId( const G4int inputValue ) { fPrimaryParticleId = inputValue; }
    void SetPrimaryParticleInitialKineticEnergy( const G4double inputValue )
      { fPrimaryParticleInitialKineticEnergy = inputValue; }
    void SetPrimaryParticleInitialTotalEnergy( const G4double inputValue )
      { fPrimaryParticleInitialTotalEnergy = inputValue; }
    void SetPrimaryParticleInitialMomentum( const G4double inputValue )
      { fPrimaryParticleInitialMomentum = inputValue; }
    void SetPrimaryParticleInitialBeta( const G4double inputValue )
      { fPrimaryParticleInitialBeta = inputValue; }
    void SetPrimaryParticleInitialGamma( const G4double inputValue )
      { fPrimaryParticleInitialGamma = inputValue; }
    void SetPrimaryParticleInitial3Momentum( const G4ThreeVector& inputValue )
      { fPrimaryParticleInitial3Momentum = inputValue; }
    void SetPrimaryParticleInitialPosition( const G4ThreeVector& inputValue )
      { fPrimaryParticleInitialPosition = inputValue; }
    void SetToleranceEPviolations( const G4double inputValue )
      { fToleranceEPviolations = inputValue; }
    void SetToleranceDeltaDecayRadius( const G4double inputValue )
      { fToleranceDeltaDecayRadius = inputValue; }
    void SetIsPreassignedDecayEnabled( const G4bool inputValue )
      { fIsPreassignedDecayEnabled = inputValue; }
    void SetIsBoostToLabEnabled( const G4bool inputValue ) { fIsBoostToLabEnabled = inputValue; }

    G4int    GetPrimaryParticleId() const { return fPrimaryParticleId; }
    G4double GetPrimaryParticleInitialKineticEnergy() const
      { return fPrimaryParticleInitialKineticEnergy; }
    G4double GetPrimaryParticleInitialTotalEnergy() const
      { return fPrimaryParticleInitialTotalEnergy; }
    G4double GetPrimaryParticleInitialMomentum() const
      { return fPrimaryParticleInitialMomentum; }
    G4double GetPrimaryParticleInitialBeta() const { return fPrimaryParticleInitialBeta; }
    G4double GetPrimaryParticleInitialGamma() const { return fPrimaryParticleInitialGamma; }
    G4ThreeVector GetPrimaryParticleInitial3Momentum() const
      { return fPrimaryParticleInitial3Momentum; }
    G4ThreeVector GetPrimaryParticleInitialPosition() const
      { return fPrimaryParticleInitialPosition; }
    G4double GetToleranceEPviolations() const { return fToleranceEPviolations; }
    G4double GetToleranceDeltaDecayRadius() const { return fToleranceDeltaDecayRadius; }
    G4bool   GetIsPreassignedDecayEnabled() const { return fIsPreassignedDecayEnabled; }
    G4bool   GetIsBoostToLabEnabled() const { return fIsBoostToLabEnabled; }

    void IncrementNumberDecays() { ++fNumDecays; }
    void IncrementNumberBadPrimaryDecays() { ++fNumBadDecays; }
    void IncrementNumberUnexpectedDecays() { ++fNumUnexpectedDecays; }
    void IncrementNumberEviolations() { ++fNumEviolations; }
    void IncrementNumberPviolations() { ++fNumPviolations; }
    void IncrementNumber_mc_truth_rPos_deltaMax_above() { ++fNum_mc_truth_rPos_deltaMax_above; }
    void IncrementNumber_underestimated_mc_truth_rPos_delta_above()
      { ++fNum_underestimated_mc_truth_rPos_delta_above; }
    void IncrementNumber_overestimated_mc_truth_rPos_delta_above()
      { ++fNum_overestimated_mc_truth_rPos_delta_above; }
    void IncrementNumberLargeUnderestimates() { ++fNumLargeUnderestimates; }
    void IncrementNumberLargeOverestimates()  { ++fNumLargeOverestimates; } 
  
    G4int GetNumberDecays() const { return fNumDecays; };
    G4int GetNumberBadDecays() const { return fNumBadDecays; }
    G4int GetNumberUnexpectedDecays() const { return fNumUnexpectedDecays; };
    G4int GetNumberEviolations() const { return fNumEviolations; };
    G4int GetNumberPviolations() const { return fNumPviolations; };
    G4int GetNumber_mc_truth_rPos_deltaMax_above() const
      { return fNum_mc_truth_rPos_deltaMax_above; }
    G4int GetNumberUnderestimated_mc_truth_rPos_delta_above() const
      { return fNum_underestimated_mc_truth_rPos_delta_above; }
    G4int GetNumberOverestimated_mc_truth_rPos_delta_above()  const
      { return fNum_overestimated_mc_truth_rPos_delta_above; }
    G4int GetNumberLargeUnderestimates() const { return fNumLargeUnderestimates; }
    G4int GetNumberLargeOverestimates()  const { return fNumLargeOverestimates; }

    void SetDecayT( const G4double inputValue );
    void SetDecayR_mc_truth( const G4double inputValue );
    void SetDecayR( const G4double inputValue );
    void SetDecayX( const G4double inputValue );
    void SetDecayY( const G4double inputValue );
    void SetDecayZ( const G4double inputValue );
    void SetDeltaDecayR( const G4double inputValue );
    void SetDeflectionAngle( const G4double inputValue );
    void SetDeltaEkin( const G4double inputValue );
    void SetDecayEkin( const G4double inputValue );
    void SetDecayPx( const G4double inputValue );
    void SetDecayPy( const G4double inputValue );
    void SetDecayPz( const G4double inputValue );
    void SetDecayEtotViolation( const G4double inputValue );
    void SetDecayPxViolation( const G4double inputValue );
    void SetDecayPyViolation( const G4double inputValue );
    void SetDecayPzViolation( const G4double inputValue );
    void SetMaxEkin_deltaMax( const G4double inputValue );
    void SetMaxEtot_deltaMax( const G4double inputValue );
    void SetMaxP_deltaMax( const G4double inputValue );
    void SetMaxPdir_deltaMax( const G4double inputValue );
    void SetMaxMass_deltaMax1( const G4double inputValue );
    void SetMaxMass_deltaMax2( const G4double inputValue );
    void SetMaxMass_deltaMax3( const G4double inputValue );
    void SetMaxBeta_deltaMax1( const G4double inputValue );
    void SetMaxBeta_deltaMax2( const G4double inputValue );
    void SetMaxGamma_deltaMax1( const G4double inputValue );
    void SetMaxGamma_deltaMax2( const G4double inputValue );
    void SetMaxGamma_deltaMax3( const G4double inputValue );
    void SetMaxT_proper_deltaMax( const G4double inputValue );
    void SetMaxT_lab_deltaMax( const G4double inputValue );
    void SetMaxMc_truth_rPos_deltaMax( const G4double inputValue );
    void SetMinUnderestimated_mc_truth_rPos_delta( const G4double inputValue );
    void SetMaxOverestimated_mc_truth_rPos_delta(  const G4double inputValue );
    void SetMinUnderestimated_rDeltaPos( const G4double inputValue );
    void SetMaxOverestimated_rDeltaPos(  const G4double inputValue );
    void SetMaxFloat_rDeltaPos_deltaMax( const G4double inputValue );

    G4double GetSumDecayT() const { return fSumDecayT; }
    G4double GetMinDecayT() const { return fMinDecayT; }
    G4double GetMaxDecayT() const { return fMaxDecayT; }
    G4double GetSumDecayR() const { return fSumDecayR; }
    G4double GetMinDecayR() const { return fMinDecayR; }
    G4double GetMaxDecayR() const { return fMaxDecayR; }
    G4double GetSumDecayX() const { return fSumDecayX; }
    G4double GetMinDecayX() const { return fMinDecayX; }
    G4double GetMaxDecayX() const { return fMaxDecayX; }
    G4double GetSumDecayY() const { return fSumDecayY; }
    G4double GetMinDecayY() const { return fMinDecayY; }
    G4double GetMaxDecayY() const { return fMaxDecayY; }
    G4double GetSumDecayZ() const { return fSumDecayZ; }
    G4double GetMinDecayZ() const { return fMinDecayZ; }
    G4double GetMaxDecayZ() const { return fMaxDecayZ; }
    G4double GetSumDeltaDecayR() const { return fSumDeltaDecayR; }
    G4double GetMinDeltaDecayR() const { return fMinDeltaDecayR; }
    G4double GetMaxDeltaDecayR() const { return fMaxDeltaDecayR; }
    G4double GetSumDeflectionAngle() const { return fSumDeflectionAngle; }
    G4double GetMinDeflectionAngle() const { return fMinDeflectionAngle; }
    G4double GetMaxDeflectionAngle() const { return fMaxDeflectionAngle; }
    G4double GetSumDeltaEkin() const { return fSumDeltaEkin; }
    G4double GetMinDeltaEkin() const { return fMinDeltaEkin; }
    G4double GetMaxDeltaEkin() const { return fMaxDeltaEkin; }
    G4double GetSumDecayEkin() const { return fSumDecayEkin; }
    G4double GetMinDecayEkin() const { return fMinDecayEkin; }
    G4double GetMaxDecayEkin() const { return fMaxDecayEkin; }
    G4double GetSumDecayPx() const { return fSumDecayPx; }
    G4double GetMinDecayPx() const { return fMinDecayPx; }
    G4double GetMaxDecayPx() const { return fMaxDecayPx; }
    G4double GetSumDecayPy() const { return fSumDecayPy; }
    G4double GetMinDecayPy() const { return fMinDecayPy; }
    G4double GetMaxDecayPy() const { return fMaxDecayPy; }
    G4double GetSumDecayPz() const { return fSumDecayPz; }
    G4double GetMinDecayPz() const { return fMinDecayPz; }
    G4double GetMaxDecayPz() const { return fMaxDecayPz; }
    G4double GetSumDecayEtotViolation() const { return fSumDecayEtotViolation; }
    G4double GetMinDecayEtotViolation() const { return fMinDecayEtotViolation; }
    G4double GetMaxDecayEtotViolation() const { return fMaxDecayEtotViolation; }
    G4double GetSumDecayPxViolation() const { return fSumDecayPxViolation; }
    G4double GetMinDecayPxViolation() const { return fMinDecayPxViolation; }
    G4double GetMaxDecayPxViolation() const { return fMaxDecayPxViolation; }
    G4double GetSumDecayPyViolation() const { return fSumDecayPyViolation; }
    G4double GetMinDecayPyViolation() const { return fMinDecayPyViolation; }
    G4double GetMaxDecayPyViolation() const { return fMaxDecayPyViolation; }
    G4double GetSumDecayPzViolation() const { return fSumDecayPzViolation; }
    G4double GetMinDecayPzViolation() const { return fMinDecayPzViolation; }
    G4double GetMaxDecayPzViolation() const { return fMaxDecayPzViolation; }
    G4double GetMaxEkin_deltaMax() const { return fMaxEkin_deltaMax; }
    G4double GetMaxEtot_deltaMax() const { return fMaxEtot_deltaMax; }
    G4double GetMaxP_deltaMax() const { return fMaxP_deltaMax; }
    G4double GetMaxPdir_deltaMax() const { return fMaxPdir_deltaMax; }
    G4double GetMaxMass_deltaMax1() const { return fMaxMass_deltaMax1; }
    G4double GetMaxMass_deltaMax2() const { return fMaxMass_deltaMax2; }
    G4double GetSumMass_deltaMax3() const { return fSumMass_deltaMax3; }
    G4double GetMaxMass_deltaMax3() const { return fMaxMass_deltaMax3; }
    G4double GetMaxBeta_deltaMax1() const { return fMaxBeta_deltaMax1; }
    G4double GetMaxBeta_deltaMax2() const { return fMaxBeta_deltaMax2; }
    G4double GetMaxGamma_deltaMax1() const { return fMaxGamma_deltaMax1; }
    G4double GetMaxGamma_deltaMax2() const { return fMaxGamma_deltaMax2; }
    G4double GetMaxGamma_deltaMax3() const { return fMaxGamma_deltaMax3; }
    G4double GetMaxT_proper_deltaMax() const { return fMaxT_proper_deltaMax; }
    G4double GetMaxT_lab_deltaMax() const { return fMaxT_lab_deltaMax; }
    G4double GetSumMc_truth_rPos_deltaMax() const { return fSumMc_truth_rPos_deltaMax; }
    G4double GetMaxMc_truth_rPos_deltaMax() const { return fMaxMc_truth_rPos_deltaMax; }
    G4double GetSumUnderestimated_mc_truth_rPos_delta() const
      { return fSumUnderestimated_mc_truth_rPos_delta; }
    G4double GetMinUnderestimated_mc_truth_rPos_delta() const
      { return fMinUnderestimated_mc_truth_rPos_delta; }
    G4double GetSumOverestimated_mc_truth_rPos_delta()  const
      { return fSumOverestimated_mc_truth_rPos_delta; }
    G4double GetMaxOverestimated_mc_truth_rPos_delta()  const
      { return fMaxOverestimated_mc_truth_rPos_delta; }
    G4double GetSumUnderestimated_rDeltaPos() const { return fSumUnderestimated_rDeltaPos; }
    G4double GetMinUnderestimated_rDeltaPos() const { return fMinUnderestimated_rDeltaPos; }
    G4double GetSumOverestimated_rDeltaPos()  const { return fSumOverestimated_rDeltaPos; }
    G4double GetMaxOverestimated_rDeltaPos()  const { return fMaxOverestimated_rDeltaPos; }
    G4double GetMaxFloat_rDeltaPos_deltaMax() const { return fMaxFloat_rDeltaPos_deltaMax; }
  
  private:
  
    G4int fNumEvents;

    G4int fPrimaryParticleId;
    G4double fPrimaryParticleInitialKineticEnergy;
    G4double fPrimaryParticleInitialTotalEnergy;
    G4double fPrimaryParticleInitialMomentum;
    G4double fPrimaryParticleInitialBeta;
    G4double fPrimaryParticleInitialGamma;
    G4ThreeVector fPrimaryParticleInitial3Momentum;
    G4ThreeVector fPrimaryParticleInitialPosition;
    G4double fToleranceEPviolations;
    G4double fToleranceDeltaDecayRadius;
    G4bool fIsPreassignedDecayEnabled;
    G4bool fIsBoostToLabEnabled;

    G4int fNumDecays;
    G4int fNumBadDecays;
    G4int fNumUnexpectedDecays;
    G4int fNumEviolations;
    G4int fNumPviolations;
    G4int fNum_mc_truth_rPos_deltaMax_above;
    G4int fNum_underestimated_mc_truth_rPos_delta_above;
    G4int fNum_overestimated_mc_truth_rPos_delta_above;
    G4int fNumLargeUnderestimates;
    G4int fNumLargeOverestimates;
  
    G4double fDecayT;
    G4double fSumDecayT;
    G4double fMinDecayT;
    G4double fMaxDecayT;
    G4double fDecayR_mc_truth;
    G4double fDecayR;
    G4double fSumDecayR;
    G4double fMinDecayR;
    G4double fMaxDecayR;
    G4double fDecayX;
    G4double fSumDecayX;
    G4double fMinDecayX;
    G4double fMaxDecayX;
    G4double fDecayY;
    G4double fSumDecayY;
    G4double fMinDecayY;
    G4double fMaxDecayY;
    G4double fDecayZ;
    G4double fSumDecayZ;
    G4double fMinDecayZ;
    G4double fMaxDecayZ;
    G4double fDeltaDecayR;
    G4double fSumDeltaDecayR;
    G4double fMinDeltaDecayR;
    G4double fMaxDeltaDecayR;
    G4double fDeflectionAngle;  // in degrees
    G4double fSumDeflectionAngle;
    G4double fMinDeflectionAngle;
    G4double fMaxDeflectionAngle;
    G4double fDeltaEkin;
    G4double fSumDeltaEkin;
    G4double fMinDeltaEkin;
    G4double fMaxDeltaEkin;
    G4double fDecayEkin;
    G4double fSumDecayEkin;
    G4double fMinDecayEkin;
    G4double fMaxDecayEkin;
    G4double fDecayPx;
    G4double fSumDecayPx;
    G4double fMinDecayPx;
    G4double fMaxDecayPx;
    G4double fDecayPy;
    G4double fSumDecayPy;
    G4double fMinDecayPy;
    G4double fMaxDecayPy;
    G4double fDecayPz;
    G4double fSumDecayPz;
    G4double fMinDecayPz;
    G4double fMaxDecayPz;
    G4double fDecayEtotViolation;
    G4double fSumDecayEtotViolation;
    G4double fMinDecayEtotViolation;
    G4double fMaxDecayEtotViolation;
    G4double fDecayPxViolation;
    G4double fSumDecayPxViolation;
    G4double fMinDecayPxViolation;
    G4double fMaxDecayPxViolation;
    G4double fDecayPyViolation;
    G4double fSumDecayPyViolation;
    G4double fMinDecayPyViolation;
    G4double fMaxDecayPyViolation;
    G4double fDecayPzViolation;
    G4double fSumDecayPzViolation;
    G4double fMinDecayPzViolation;
    G4double fMaxDecayPzViolation;
    G4double fMaxEkin_deltaMax;
    G4double fMaxEtot_deltaMax;
    G4double fMaxP_deltaMax;
    G4double fMaxPdir_deltaMax;
    G4double fMaxMass_deltaMax1;
    G4double fMaxMass_deltaMax2;
    G4double fSumMass_deltaMax3;
    G4double fMaxMass_deltaMax3;
    G4double fMaxBeta_deltaMax1;
    G4double fMaxBeta_deltaMax2;
    G4double fMaxGamma_deltaMax1;
    G4double fMaxGamma_deltaMax2;
    G4double fMaxGamma_deltaMax3;
    G4double fMaxT_proper_deltaMax;
    G4double fMaxT_lab_deltaMax;
    G4double fSumMc_truth_rPos_deltaMax;
    G4double fMaxMc_truth_rPos_deltaMax;
    G4double fSumUnderestimated_mc_truth_rPos_delta;
    G4double fMinUnderestimated_mc_truth_rPos_delta;
    G4double fSumOverestimated_mc_truth_rPos_delta;
    G4double fMaxOverestimated_mc_truth_rPos_delta;
    G4double fSumUnderestimated_rDeltaPos;
    G4double fMinUnderestimated_rDeltaPos;
    G4double fSumOverestimated_rDeltaPos;  
    G4double fMaxOverestimated_rDeltaPos;
    G4double fMaxFloat_rDeltaPos_deltaMax;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
