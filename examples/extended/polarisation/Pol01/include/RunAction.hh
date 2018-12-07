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
/// \file polarisation/Pol01/include/RunAction.hh
/// \brief Definition of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "ProcessesCount.hh"
#include "globals.hh"
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class PrimaryGeneratorAction;
class G4Run;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  class ParticleStatistics {
  public:
    ParticleStatistics();
    ~ParticleStatistics();
    void EventFinished();
    void FillData(G4double kinEnergy, G4double costheta,
                  G4double longitudinalPolarization);
    void PrintResults(G4int totalNumberOfEvents);
    void Clear();
  private:
    G4int fCurrentNumber;
    G4int fTotalNumber, fTotalNumber2;
    G4double fSumEnergy, fSumEnergy2;
    G4double fSumPolarization, fSumPolarization2;
    G4double fSumCosTheta, fSumCosTheta2;
  };

public:

  RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);

  void CountProcesses(G4String);

  void FillData(const G4ParticleDefinition* particle,
                G4double kinEnergy, G4double costheta, G4double phi,
                G4double longitudinalPolarization);

  void EventFinished();
                                     
private:

  void BookHisto();
  void SaveHisto(G4int nevents);

  const G4ParticleDefinition* fGamma;
  const G4ParticleDefinition* fElectron;
  const G4ParticleDefinition* fPositron;

  DetectorConstruction*   fDetector;
  PrimaryGeneratorAction* fPrimary;
  ProcessesCount*         fProcCounter;

  G4AnalysisManager*      fAnalysisManager;
  
  G4int fTotalEventCount;

  ParticleStatistics fPhotonStats;
  ParticleStatistics fElectronStats;
  ParticleStatistics fPositronStats;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
