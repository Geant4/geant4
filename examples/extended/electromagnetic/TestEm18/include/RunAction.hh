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
/// \file electromagnetic/TestEm18/include/RunAction.hh
/// \brief Definition of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class G4ParticleDefinition;
class G4Material;

class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void CountProcesses(G4String procName);

    void TrackLength (G4double step);

    void EnergyDeposited (G4double edepPrim, G4double edepSecond);

    void EnergyTransferedByProcess (G4String procName, G4double energy);

    void EnergyTransfered (G4double energy);

    void TotalEnergyLost (G4double energy);

    void EnergyBalance (G4double energy);

    void TotalEnergyDeposit (G4double energy);

    void EnergySpectrumOfSecondaries (G4String particleName, G4double ekin);

  public:
    G4double GetEnergyFromRestrictedRange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);
                       
    G4double GetEnergyFromCSDARange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);

private:
  struct MinMaxData {
   MinMaxData()
     : fCount(0), fVsum(0.), fVmin(0.), fVmax(0.) {}
   MinMaxData(G4int count, G4double vsum, G4double vmin, G4double vmax)
     : fCount(count), fVsum(vsum), fVmin(vmin), fVmax(vmax) {}
   G4int     fCount;
   G4double  fVsum;
   G4double  fVmin;
   G4double  fVmax;
  };
  
  private:

    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    HistoManager*           fHistoManager;

    std::map<G4String,G4int>  fProcCounter;

    G4long   fNbSteps;
    G4double fTrackLength, fStepMin, fStepMax;

    G4double fEdepPrimary, fEdepPrimMin, fEdepPrimMax;
    std::map<G4String,MinMaxData> fEtransfByProcess;
    G4double fEnergyTransfered, fEtransfMin, fEtransfMax;
    G4double fEnergyLost, fElostMin, fElostMax;
    G4double fEnergyBalance, fEbalMin, fEbalMax;

    G4double fEdepSecondary, fEdepSecMin, fEdepSecMax;
    G4double fEdepTotal, fEdepTotMin, fEdepTotMax;

    std::map<G4String,MinMaxData> fEkinOfSecondaries;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

