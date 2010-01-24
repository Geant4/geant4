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
// $Id: RunAction.hh,v 1.20 2010-01-24 17:25:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  RunAction(DetectorConstruction*, PrimaryGeneratorAction*, HistoManager*);
  ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);

  void fillPerEvent(G4int,G4double,G4double);
    
  void sumEnergyFlow(G4int plane, G4double Eflow)
                                            {EnergyFlow[plane]  += Eflow;};
  void sumLateralEleak(G4int cell, G4double Eflow)
                                            {lateralEleak[cell] += Eflow;};
    
  void PrintDedxTables();

  void AddSecondaryTrack(const G4Track*);
    
  // Acceptance parameters
  void     SetEdepAndRMS(G4int, G4double, G4double, G4double);
  G4double GetAverageEdep(G4int i) const    {return edeptrue[i];};
  G4double GetRMSEdep(G4int i) const        {return rmstrue[i];};
  G4double GetLimitEdep(G4int i) const      {return limittrue[i];};
  void     SetApplyLimit(G4bool val)        {applyLimit = val;};

private:
  
  DetectorConstruction*   Detector;
  PrimaryGeneratorAction* Primary;    
  RunActionMessenger*     runMessenger;
  HistoManager*           histoManager;

  G4double sumEAbs [MaxAbsor], sum2EAbs [MaxAbsor]; 
  G4double sumLAbs [MaxAbsor], sum2LAbs [MaxAbsor];
    
  std::vector<G4double> EnergyFlow;
  std::vector<G4double> lateralEleak;
  std::vector<G4double> energyDeposit[MaxAbsor];
    
  G4double edeptrue [MaxAbsor];
  G4double rmstrue  [MaxAbsor];
  G4double limittrue[MaxAbsor];

  G4int  n_gamma;
  G4int  n_elec;
  G4int  n_pos;

  G4bool applyLimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

