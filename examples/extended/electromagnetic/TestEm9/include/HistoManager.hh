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
/// \file electromagnetic/TestEm9/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//

#ifndef HistoManager_h
#define HistoManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to hold Emc parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 27/09/00
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DataVector.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;

class HistoManager
{

public:
  // With description

  static HistoManager* GetPointer();

private:

  HistoManager();

public: // Without description

  ~HistoManager();

  void BookHisto();

  void BeginOfRun();
  void EndOfRun(G4int runID);

  void BeginOfEvent();
  void EndOfEvent();

  void ScoreNewTrack(const G4Track* aTrack);
  void AddEnergy(G4double edep, G4int idx, G4int copyNo);

  void AddDeltaElectron(const G4DynamicParticle*);
  void AddPhoton(const G4DynamicParticle*);

  inline void ResetTrackLength() {fTrackLength = 0.0, fTrackAbs = true;};
  inline void AddPositron(const G4DynamicParticle*) {++fPosit;};
  inline void SetVerbose(G4int val) {fVerbose = val;};
  inline G4int GetVerbose() const {return fVerbose;};
  inline void SetHistoNumber(G4int val) {fNHisto = val;};

  inline void SetFirstEventToDebug(G4int val) {fEvt1 = val;};
  inline G4int FirstEventToDebug() const {return fEvt1;};
  inline void SetLastEventToDebug(G4int val) {fEvt2 = val;};
  inline G4int LastEventToDebug() const {return fEvt2;};

  inline void SetMaxEnergy(G4double val) {fMaxEnergy = val;};
  inline G4double  GetMaxEnergy() const {return fMaxEnergy;};
  inline void AddStep() {fStep += 1.0;};

  // Acceptance parameters
  inline void SetBeamEnergy(G4double val) {fBeamEnergy = val;};

  void SetEdepAndRMS(G4int, G4ThreeVector);

private:

  // MEMBERS
  static HistoManager* fManager;

  const G4ParticleDefinition* fGamma;
  const G4ParticleDefinition* fElectron;
  const G4ParticleDefinition* fPositron;

  G4int fNHisto;
  G4int fVerbose;
  G4int fEvt1;
  G4int fEvt2;

  G4double fBeamEnergy;
  G4double fMaxEnergy;
  G4double fMaxEnergyAbs;

  G4double fTrackLength;
  G4double fStep;
//  G4double fStepTarget;
  G4bool fTrackAbs;        // Track is in absorber
  G4int fEvt;
  G4int fElec;
  G4int fPosit;
  G4int fGam;
//  G4int fGamph;
//  G4int fGamTar;
  G4int fLowe;
  G4int fBinsE, fBinsEA, fBinsED;

  G4double fEabs1, fEabs2, fEabs3, fEabs4;
  G4double     fE[25];
  G4DataVector fEvertex;
  G4DataVector fNvertex;
  G4DataVector fBrem;
  G4DataVector fPhot;
  G4DataVector fComp;
  G4DataVector fConv;

  G4double  fEdeptrue[3];
  G4double  fRmstrue[3];
  G4double  fLimittrue[3];
  G4double  fEdep[6];
  G4double  fErms[6];
  G4double  fEdeptr[6];
  G4double  fErmstr[6];
  G4int     fStat[6];
  G4int     fNmax;

  Histo*    fHisto;
  


};

#endif
