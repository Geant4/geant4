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
// This example is provided by the Geant4-DNA collaboration
// dnadamage3 example is derived from the chem6 example
// chem6 example authors: W. G. Shin and S. Incerti (CENBG, France)
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file ScoreStrandBreaks.hh
/// \brief Definition of the ScoreStrandBreaks class

#ifndef DNADAMAGE2_ScoreSB_h
#define DNADAMAGE2_ScoreSB_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImessenger.hh"
#include "G4DNAIRT.hh"
#include "DetectorConstruction.hh"
#include "MoleculeInserter.hh"

class G4VAnalysisManager;
class G4MolecularConfiguration;

typedef std::map<G4int,std::map<G4int,std::map<G4int,
                       std::map<G4int, G4double>>>> EnergyDepositMap;
typedef std::map<G4int,std::vector<std::vector<G4int>>> DamageMap;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ScoreStrandBreaks : public G4VPrimitiveScorer,
                          public G4UImessenger
{
public:
  ScoreStrandBreaks(G4String name,DetectorConstruction*,G4double*);
  ~ScoreStrandBreaks() override;

  void WriteWithAnalysisManager(G4VAnalysisManager*);
  DamageMap GetDamageMap() { return fDirectDamageMap; }
  void ClearDamageMap() { fDirectDamageMap.clear(); }

  void Initialize(G4HCofThisEvent*) override;
  void EndOfEvent(G4HCofThisEvent*) override;
  void Clear();
  void DrawAll() override;
  void PrintAll() override;
  void AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* );
  void OutputAndClear(G4double, G4double);
  void SetNewValue(G4UIcommand*, G4String) override;

  void ASCII(G4double, G4double);
  
private:
  G4int    fnbOfEvents = 0;
  G4bool   fDNAInserted = false;
  G4double fBreakEnergy = 17.5 * eV;
  G4double fEnergyDeposit = 0;
  G4double* fRadius = nullptr;

  G4String fOutputName = "DirectDamageInfo";
  G4String fOutputType = "ascii";
  DamageMap  fDirectDamageMap;
  DamageMap  fIndirectDamageMap;
  EnergyDepositMap fEnergyDepositMap;

  G4UIcmdWithAString* fpOutputFileUI = nullptr;
  G4UIcmdWithAString* fpOutputTypeUI = nullptr;
  G4UIcmdWithADoubleAndUnit* fpBreakEnergyUI = nullptr;

  DetectorConstruction* fpDetector = nullptr;
  std::map<G4int,G4double> fDoseArray;

  MoleculeInserter* fpGun = nullptr;

protected:
  G4bool ProcessHits(G4Step*,G4TouchableHistory*) override;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
