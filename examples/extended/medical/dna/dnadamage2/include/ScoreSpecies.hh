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
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file ScoreSpecies.hh
/// \brief Definition of the ScoreSpecies class
// Description:
//   This is a primitive scorer class for scoring the radiolitic species
// produced after irradiation in a water volume
//
// Created: 2015-10-27  by M. Karamitros,
// modified: 2016-03-16 by P. Piersimoni
// modified: 2022-01-2022 by J. Naoki D. Kondo

#ifndef DNADAMAGE2_ScoreSpecies_h
#define DNADAMAGE2_ScoreSpecies_h 1

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

class ScoreSpecies : public G4VPrimitiveScorer,
                     public G4UImessenger
{
public:
  ScoreSpecies(G4String name, G4int depth=0);
  ~ScoreSpecies() override;
  
  /** Add a time at which the number of species should be recorded.
      Default times are set up to 1 microsecond.*/
  inline void AddTimeToRecord(double time)
  {
    fTimeToRecord.insert(time);
  }
  
  /**  Remove all times to record, must be reset by user.*/
  inline void ClearTimeToRecord()
  {
    fTimeToRecord.clear();
  }

  /** Get number of recorded events*/
  inline int GetNumberOfRecordedEvents() const
  {
    return fNEvent;
  }

  /** Write results to whatever chosen file format*/
  void WriteWithAnalysisManager(G4VAnalysisManager*);
  
  struct SpeciesInfo
  {
    SpeciesInfo() {}
    SpeciesInfo(const SpeciesInfo& right) // Species A(B);
    {
      fNumber  = right.fNumber;
      fNumber2 = right.fNumber2;
      fG = right.fG;
      fG2 = right.fG2;
    }
    SpeciesInfo& operator=(const SpeciesInfo& right) // A = B
    {
      if(&right == this) return *this;
      fNumber2 = right.fNumber2;
      fNumber  = right.fNumber;
      fG  = right.fG;
      fG2 = right.fG2;
      return *this;
    }
    int fNumber   = 0;
    int fNumber2  = 0;
    double fG  = 0;
    double fG2 = 0;
  };
  
  
private:
  typedef const G4MolecularConfiguration Species;
  typedef std::map<Species*, SpeciesInfo>  InnerSpeciesMap;
  typedef std::map<double, InnerSpeciesMap> SpeciesMap;
  SpeciesMap fSpeciesInfoPerTime;

  std::set<G4double> fTimeToRecord;
  
  int fNEvent = 0; // number of processed events
  double fEdep = 0; // total energy deposition
  G4String fOutputType; // output type 

protected:
  G4bool ProcessHits(G4Step*,G4TouchableHistory*) override;

public:
  void Initialize(G4HCofThisEvent*) override;
  void EndOfEvent(G4HCofThisEvent*) override;
  void DrawAll() override;
  void PrintAll() override;
  /** Method used in multithreading mode in order to merge 
      the results*/
  void AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* );
  void OutputAndClear();
  void SetNewValue(G4UIcommand*, G4String) override;
  void OutputToASCII();
  
  SpeciesMap GetSpeciesInfo() {return fSpeciesInfoPerTime;}

private:
  G4int fHCID;
  G4THitsMap<G4double>* fEvtMap = nullptr;

  G4int fRunID = 0;
  G4UIdirectory* fSpeciesdir = nullptr;
  G4UIcmdWithAnInteger* fTimeBincmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fAddTimeToRecordcmd = nullptr;

  G4UIcmdWithAString* fOutputTypeUI = nullptr;
  G4UIcmdWithAString* fOutputFileUI = nullptr;

  G4String fOutputFile = "Species_Info";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
