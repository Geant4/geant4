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
/// \file ScoreSpecies.hh
/// \brief Definition of the ScoreSpecies class

// This example is provided by the Geant4-DNA collaboration
// chem6 example is derived from chem4 and chem5 examples
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: W. G. Shin and S. Incerti (CENBG, France)

#ifndef CHEM6_ScoreSpecies_h
#define CHEM6_ScoreSpecies_h 1

#include "G4THitsMap.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UImessenger.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4UIcmdWithAString.hh"
#include <set>
#include <memory>


class G4VAnalysisManager;
class AnalysisMessenger;
class G4MolecularConfiguration;

/** \file ScoreSpecies.hh*/

// Description:
//   This is a primitive scorer class for scoring the radiolitic species
// produced after irradiation in a water volume
// from chem4 example

class ScoreSpecies final : public G4VPrimitiveScorer, public G4UImessenger
{
    struct SpeciesInfo
    {
      SpeciesInfo()
      {
        fNumber = 0;
        fG = 0.;
        fG2 = 0.;
      }
      SpeciesInfo(const SpeciesInfo& right)  // Species A(B);
      {
        fNumber = right.fNumber;
        fG = right.fG;
        fG2 = right.fG2;
      }
      SpeciesInfo& operator=(const SpeciesInfo& right)  // A = B
      {
        if (&right == this) return *this;
        fNumber = right.fNumber;
        fG = right.fG;
        fG2 = right.fG2;
        return *this;
      }
      G4int fNumber;
      G4double fG;
      G4double fG2;
    };
    using Species = const G4MolecularConfiguration;
    using InnerSpeciesMap = std::map<Species*, SpeciesInfo>;
    using SpeciesMap =  std::map<G4double, InnerSpeciesMap>;

  public:
    explicit ScoreSpecies(const G4String &name, G4int depth = 0);

    ~ScoreSpecies() override = default;

      /** Add a time at which the number of species should be recorded.
          Default times are set up to 1 microsecond.*/
    void AddTimeToRecord(const G4double time) { fTimeToRecord.insert(time); }

      /**  Remove all times to record, must be reset by user.*/
    void ClearTimeToRecord() { fTimeToRecord.clear(); }

      /** Get number of recorded events*/
    G4int GetNumberOfRecordedEvents() const { return fNEvent; }

      /** Write results to whatever chosen file format*/
    void WriteWithAnalysisManager(G4VAnalysisManager*);

    void Initialize(G4HCofThisEvent*) override;
    void EndOfEvent(G4HCofThisEvent*) override;
    void DrawAll() override{};
    void PrintAll() override;
    /** Method used in multithreading mode in order to merge
        the results*/
    void AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer*);
    void OutputAndClear();
    void SetNewValue(G4UIcommand*, G4String) override;

    void SetFileName(const G4String& name) { fFileName = name; };
    SpeciesMap GetSpeciesInfo() { return fSpeciesInfoPerTime; }

  protected:
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
  private:
    SpeciesMap fSpeciesInfoPerTime;
    std::set<G4double> fTimeToRecord;
    G4int fNEvent = 0;  // number of processed events
    G4double fEdep = 0;  // total energy deposition
    G4String fOutputType = "root";  // output type
    G4int fHCID = -1;
    G4THitsMap<G4double>* fEvtMap = nullptr;
    G4int fRunID = 0;
    G4String fFileName = "Species";
    std::unique_ptr<G4UIdirectory> fSpeciesdir;
    std::unique_ptr<G4UIcmdWithAnInteger> fTimeBincmd;
    std::unique_ptr<G4UIcmdWithADoubleAndUnit> fAddTimeToRecordcmd;
    std::unique_ptr<G4UIcmdWithAString> fNameCmd;
};
#endif
