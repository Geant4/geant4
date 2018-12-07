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
/// \file electromagnetic/TestEm8/include/Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4DataVector.hh"
#include "G4StatDouble.hh"

#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Step;
class G4ElectronIonPair;
class TestParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
public:

  Run();
  virtual ~Run();

  virtual void Merge(const G4Run*);

  void BeginOfRun();
  void EndOfRun();

  void BeginOfEvent();
  void EndOfEvent();

  void AddEnergy(G4double edep, const G4Step*);

  inline void SetVerbose(G4int value);

  inline G4int GetVerbose() const;

  inline G4double GetTotStepGas() const;

  inline G4double GetTotCluster() const;

  inline G4double GetMeanCluster() const;

  inline const G4StatDouble* GetStat() const;

private:

  G4int fVerbose;
  G4int fNbins;
  G4double fStepGas;
  G4double fMaxEnergy;
  G4double fCluster;
  G4double fTotStepGas;
  G4double fTotCluster;
  G4double fMeanCluster;
  G4double fFactorALICE;
  G4double fWidthALICE;
  G4double fEvt;

  G4double fTotEdep;
  G4StatDouble fEdep;
  G4double fOverflow;
  G4DataVector fEgas;

  G4ElectronIonPair* fElIonPair;
  TestParameters* fParam;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void Run::SetVerbose(G4int value)
{
  fVerbose = value;
}

inline G4int Run::GetVerbose() const
{
  return fVerbose;
}

inline G4double Run::GetTotStepGas() const
{
  return fTotStepGas;
}

inline G4double Run::GetTotCluster() const
{
  return fTotCluster;
} 

inline G4double Run::GetMeanCluster() const
{
  return fMeanCluster;
}

inline const G4StatDouble* Run::GetStat() const
{
  return &fEdep;
}

#endif


