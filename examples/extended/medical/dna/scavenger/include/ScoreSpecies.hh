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
/// \file scavenger/include/ScoreSpecies.hh
/// \brief Definition of the scavenger::ScoreSpecies class

#ifndef SCAVENGER_ScoreSpecies_h
#define SCAVENGER_ScoreSpecies_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include <set>
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImessenger.hh"

class G4VAnalysisManager;

class G4MolecularConfiguration;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace scavenger
{

/// \brief Primitive scorer class for scoring the radiolytic species
/// produced after irradiation in a water volume
///
/// This is a primitive scorer class for molecular species.
/// The number of species is recorded for all times (predetermined or
/// user chosen). It also scores the energy deposition in order to compute the
/// radiochemical yields.

class ScoreSpecies
  : public G4VPrimitiveScorer
  , public G4UImessenger
{
 public:
  explicit ScoreSpecies(const G4String& name, const G4int& depth = 0);

  ~ScoreSpecies() override = default;

  /** Add a time at which the number of species should be recorded.
   Default times are set up to 1 microsecond.*/
  inline void AddTimeToRecord(const G4double& time)
  {
    fTimeToRecord.insert(time);
  }

  /**  Remove all times to record, must be reset by user.*/
  inline void ClearTimeToRecord() { fTimeToRecord.clear(); }

  /** Get number of recorded events*/
  [[nodiscard]] inline G4int GetNumberOfRecordedEvents() const
  {
    return fNEvent;
  }

  void Initialize(G4HCofThisEvent*) override;

  void EndOfEvent(G4HCofThisEvent*) override;

  void DrawAll() override{};

  void PrintAll() override;

  /** Method used in multithreading mode in order to merge
   the results*/
  void AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer*);

  void OutputAndClear();

  void SetNewValue(G4UIcommand*, G4String) override;

  /** Write results to whatever chosen file format*/
  void WriteWithAnalysisManager(G4VAnalysisManager*);

  struct SpeciesInfo
  {
    SpeciesInfo() = default;
    SpeciesInfo(const SpeciesInfo& right) = default;
    SpeciesInfo& operator=(const SpeciesInfo& right) = default;

    G4int fNumber = 0;
    G4double fG   = 0.;
    G4double fG2  = 0;
  };

 protected :
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

 private:
  typedef const G4MolecularConfiguration Species;
  typedef std::map<Species*, SpeciesInfo> InnerSpeciesMap;
  typedef std::map<G4double, InnerSpeciesMap> SpeciesMap;
  SpeciesMap fSpeciesInfoPerTime;
  std::set<G4double> fTimeToRecord;
  G4int fNEvent = 0;   // number of processed events
  G4double fEdep = 0;  // total energy deposition
  G4String fOutputType = "root";
  G4int fHCID = -1;
  G4THitsMap<G4double>* fEvtMap = nullptr;
  std::unique_ptr<G4UIdirectory> fpSpeciesdir;
  std::unique_ptr<G4UIcmdWithAnInteger> fpTimeBincmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fpAddTimeToRecordcmd;
  std::unique_ptr<G4UIcmdWithAString> fpSetResultsFileNameCmd;
  G4String fRootFileName = "scorer.root";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

}