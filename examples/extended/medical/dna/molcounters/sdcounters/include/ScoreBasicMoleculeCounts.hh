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
/// \file ScoreBasicMoleculeCounts.hh
/// \brief Definition of the ScoreBasicMoleculeCounts class

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025)
//

#ifndef ScoreBasicMoleculeCounts_hh
#define ScoreBasicMoleculeCounts_hh 1

#include "G4MoleculeCounter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UImessenger.hh"
#include "G4VPrimitiveScorer.hh"

#include <set>

class G4HCofThisEvent;
class G4MolecularConfiguration;
class G4VAnalysisManager;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ScoreBasicMoleculeCounts : public G4VPrimitiveScorer, public G4UImessenger
{
  public:
    ScoreBasicMoleculeCounts(G4String name, G4int depth = 0, G4String = "BasicCounter");
    ~ScoreBasicMoleculeCounts() override;

    void EndOfEvent(G4HCofThisEvent*) override;
    void clear() override;
    void AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer*);
    void OutputAndClear();
    void SetNewValue(G4UIcommand*, G4String) override;
    const std::map<G4String, std::map<G4double, G4double>>& GetMoleculeCountMap() const
    {
      return fMoleculeCountPerIndexPerTime;
    }
    inline void AddTimeToRecord(G4double);
    inline void ClearTimesToRecord();
    inline G4int GetNumberOfRecordedEvents() const;
    void WriteWithAnalysisManager();
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

  private:
    const char* ScorerCommand(G4String command);
    G4String fMoleculeCounterName = "";
    G4int fRunID = 0;
    G4int fNbOfScoredEvents = 0;
    G4UIcmdWithAnInteger* fTimeBincmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fAddTimeToRecordcmd = nullptr;
    //G4bool fAppendAllCopyNumbers{false};//not used
    std::set<G4double> fTimesToRecord;
    std::map<G4String, std::map<G4double, G4double>> fMoleculeCountPerIndexPerTime;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void ScoreBasicMoleculeCounts::AddTimeToRecord(G4double time)
{
  fTimesToRecord.insert(time);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void ScoreBasicMoleculeCounts::ClearTimesToRecord()
{
  fTimesToRecord.clear();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline G4int ScoreBasicMoleculeCounts::GetNumberOfRecordedEvents() const
{
  return fNbOfScoredEvents;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
