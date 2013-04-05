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
/// \file electromagnetic/TestEm3/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "DetectorConstruction.hh"

#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction;
class RunActionMessenger;
class HistoManager;
class G4Track;
class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
public:

  RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
 ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);

  void FillPerEvent(G4int,G4double,G4double);
    
  void SumEnergyFlow(G4int plane, G4double Eflow)
                                            {fEnergyFlow[plane]  += Eflow;};
  void SumLateralEleak(G4int cell, G4double Eflow)
                                            {fLateralEleak[cell] += Eflow;};

  void AddChargedStep() { fChargedStep += 1.0; }
  void AddNeutralStep() { fNeutralStep += 1.0; }
    
  void PrintDedxTables();

  void AddSecondaryTrack(const G4Track*);
    
  // Acceptance parameters
  void     SetEdepAndRMS(G4int, G4double, G4double, G4double);
  G4double GetAverageEdep(G4int i) const    {return fEdeptrue[i];};
  G4double GetRMSEdep(G4int i) const        {return fRmstrue[i];};
  G4double GetLimitEdep(G4int i) const      {return fLimittrue[i];};
  void     SetApplyLimit(G4bool val)        {fApplyLimit = val;};

private:
  
  DetectorConstruction*   fDetector;
  PrimaryGeneratorAction* fPrimary;    
  RunActionMessenger*     fRunMessenger;
  HistoManager*           fHistoManager;

  G4double fSumEAbs [MaxAbsor], fSum2EAbs [MaxAbsor]; 
  G4double fSumLAbs [MaxAbsor], fSum2LAbs [MaxAbsor];
    
  std::vector<G4double> fEnergyFlow;
  std::vector<G4double> fLateralEleak;
  std::vector<G4double> fEnergyDeposit[MaxAbsor];
    
  G4double fEdeptrue [MaxAbsor];
  G4double fRmstrue  [MaxAbsor];
  G4double fLimittrue[MaxAbsor];

  G4double fChargedStep;
  G4double fNeutralStep;

  G4int  fN_gamma;
  G4int  fN_elec;
  G4int  fN_pos;

  G4bool fApplyLimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

