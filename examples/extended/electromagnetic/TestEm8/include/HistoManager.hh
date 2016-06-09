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
/// \file electromagnetic/TestEm8/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to make analysis and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

#ifndef HistoManager_h
#define HistoManager_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4DataVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;
class G4Step;
class G4ElectronIonPair;

class HistoManager
{

public:
  // With description

  static HistoManager* GetPointer();

private:

  HistoManager();

public: // Without description

  ~HistoManager();

  void bookHisto();

  void BeginOfRun();
  void EndOfRun();

  void BeginOfEvent();
  void EndOfEvent();

  void AddEnergy(G4double edep, G4Step*);

  inline void SetMaxEnergy(G4double value);

  inline void SetNumberBins(G4int value);

  inline void SetNumberBinsCluster(G4int value);

  inline void SetVerbose(G4int value);

  inline G4int GetVerbose() const;

private:

  // MEMBERS
  static HistoManager* fManager;

  G4int fNHisto;
  G4int fVerbose;

  G4double fMaxEnergy;
  G4double fStepGas;
  G4double fCluster;
  G4double fTotStepGas;
  G4double fTotCluster;
  G4double fEvt;

  G4int fBinsE; 
  G4int fBinsCluster;

  G4bool fHistoBooked;

  G4double fTotEdep;
  G4double fOverflow;
  G4DataVector fEgas;

  Histo*   fHisto;
  G4ElectronIonPair* fElIonPair;
};

inline void HistoManager::SetMaxEnergy(G4double value)
{
  fMaxEnergy = value;
}

inline void HistoManager::SetNumberBins(G4int value)
{
  fBinsE = value;
}

inline void HistoManager::SetNumberBinsCluster(G4int value)
{
  fBinsCluster = value;
}

inline void HistoManager::SetVerbose(G4int value)
{
  fVerbose = value;
}

inline G4int HistoManager::GetVerbose() const
{
  return fVerbose;
}

#endif
