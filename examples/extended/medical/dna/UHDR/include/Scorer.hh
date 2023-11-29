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
#ifndef scorer_h
#define scorer_h

#include "G4DNAMesh.hh"
#include "G4THitsMap.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"
#include "G4VPrimitiveScorer.hh"
#include <memory>
#include <set>

class G4DNAEventScheduler;

class G4VAnalysisManager;

class G4MolecularConfiguration;

class G4VChemistryWorld;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

struct Dose : public G4UImessenger {
  Dose();

  ~Dose() override = default;

  void SetNewValue(G4UIcommand *, G4String) final;

  std::unique_ptr<G4UIdirectory> fpDoseDir;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fpAddDoseCutOff;
  G4double fDosesCutOff = 0;
  G4double fCumulatedDose = 0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

struct Gvalues : public G4UImessenger {
  Gvalues();

  ~Gvalues() override = default;

  void SetNewValue(G4UIcommand *, G4String) final;

  G4int fHCID = -1;
  G4double fTimeLimit;
  std::unique_ptr<G4UIdirectory> fSpeciesdir;
  std::unique_ptr<G4UIcmdWithAnInteger> fTimeBincmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fAddTimeToRecordcmd;

  G4int fNEvent = 0;
  G4double fEdep = 0;

  inline void AddTimeToRecord(double time) { fTimeToRecord.insert(time); }

  void WriteWithAnalysisManager(G4VAnalysisManager *, const std::string &out);

  struct SpeciesInfo {
    SpeciesInfo() = default;

    SpeciesInfo(const SpeciesInfo &right) = default;

    SpeciesInfo &operator=(const SpeciesInfo &right) = default;

    int64_t fNumber = 0;
    G4double fG = 0;
    G4double fG2 = 0;
  };

  inline auto GetNumberOfRecordedEvents() const { return fNEvent; }

  using Species = const G4MolecularConfiguration;
  using InnerSpeciesMap = std::map<Species *, SpeciesInfo>;
  using SpeciesMap = std::map<double, InnerSpeciesMap>;
  SpeciesMap fSpeciesInfoPerTime;
  std::set<G4double> fTimeToRecord;
  G4int fRunID = 1;

  inline auto GetSpeciesInfo() const { return fSpeciesInfoPerTime; }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

template<typename TR>
class Scorer : public G4VPrimitiveScorer {
public: // with description
  Scorer();

  ~Scorer() override = default;

  void Initialize(G4HCofThisEvent *) override;

  void EndOfEvent(G4HCofThisEvent *) override;

  void clear() override;

  G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;

  G4THitsMap<G4double> *GetEventMap() const;

  void AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer *);

  void OutputAndClear(const std::string &out);

  TR *GetpScorer();

  void SetChemistryWorld(G4VChemistryWorld *);

  G4VChemistryWorld *GetChemistryWorld() const;

  void SaveScavengerChange();

  void SaveMoleculeCounter();

  void SetEventScheduler(G4DNAEventScheduler *pEventScheduler) {
    fpEventScheduler = pEventScheduler;
  }

private:
  std::unique_ptr<TR> fpScorer;
  G4int fHCID = -1;
  G4THitsMap<G4double> *fpEvtMap = nullptr;
  G4VChemistryWorld *fpChemistryWorld = nullptr;
  G4DNAEventScheduler *fpEventScheduler = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

template<typename TR>
Scorer<TR>::Scorer()
    : G4VPrimitiveScorer(typeid(TR).name()), fpScorer(new TR){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

template<typename TR>
TR *Scorer<TR>::GetpScorer() { return fpScorer.get(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

template<typename TR>
G4THitsMap<G4double> *Scorer<TR>::GetEventMap() const {
  return fpEvtMap;
}

#endif