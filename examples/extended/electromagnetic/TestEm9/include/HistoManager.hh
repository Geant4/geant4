//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: HistoManager.hh,v 1.3 2003-10-31 12:08:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
#include "Histo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class HistoManager
{

public:
  // With description

  static HistoManager* GetPointer();

private:

  HistoManager();

public: // Without description

  ~HistoManager();

  void BeginOfEvent();
  void EndOfEvent();

  void BeginOfRun();
  void EndOfRun();

  void bookHisto();

  void SaveToTuple(const G4String&, G4double);
  void SaveToTuple(const G4String&, G4double, G4double);
  void SaveEvent();
  G4double GetTrackLength() const {return trackLength;};
  void ResetTrackLength() {trackLength = 0.0, trackAbs = true;};
  void SetTrackOutAbsorber() {trackAbs = false;};
  G4bool GetTrackInAbsorber() const {return trackAbs;};
  void AddTrackLength(G4double x)   {trackLength += x;};
  void AddEndPoint(G4double);
  void AddEnergy(G4double, G4double);
  void AddDeltaElectron(const G4DynamicParticle*);
  void AddPhoton(const G4DynamicParticle*);
  void AddPositron(const G4DynamicParticle*) {n_posit++;};
  void SetVerbose(G4int val) {verbose = val;};
  G4int GetVerbose() const {return verbose;};
  void SetHistoNumber(G4int val) {nHisto = val;};
  void SetNtuple(G4bool val) {nTuple = val;};

  void SetFirstEventToDebug(G4int val) {nEvt1 = val;};
  G4int FirstEventToDebug() const {return nEvt1;};
  void SetLastEventToDebug(G4int val) {nEvt2 = val;};
  G4int LastEventToDebug() const {return nEvt2;};

  void SetMaxEnergy(G4double val) {maxEnergy = val;};
  G4double  GetMaxEnergy() const {return maxEnergy;};
  void SetThresholdEnergy(G4double val) {thKinE = val;};
  void SetThresholdZ(G4double val) {thPosZ = val;};
  void AddStep() {n_step++;};

  void AddEnergy(G4double edep, G4int idx, G4int copyNo);
  void ScoreNewTrack(const G4Track* aTrack);

private:

  // MEMBERS
  static HistoManager* fManager;

  G4int nHisto;
  G4int verbose;
  G4int nEvt1;
  G4int nEvt2;

  G4double beamEnergy;
  G4double maxEnergy;
  G4double maxEnergyAbs;
  G4double thKinE;
  G4double thPosZ;

  G4double trackLength;
  G4bool trackAbs;        // Track is in absorber
  G4int n_evt;
  G4int n_elec;
  G4int n_posit;
  G4int n_gam;
  G4int n_step;
  G4int n_gamph;
  G4int n_gam_tar;
  G4int n_step_target;
  G4int tCounter;
  G4int nBinsE, nBinsEA, nBinsED;
  G4bool nTuple;

  G4double Eabs1, Eabs2, Eabs3, Eabs4;
  G4double E[25];
  G4DataVector Evertex;
  G4DataVector Nvertex;

  std::vector<G4int> tScore;
  std::vector<G4int> tType;

  Histo  histo;
};

#endif
