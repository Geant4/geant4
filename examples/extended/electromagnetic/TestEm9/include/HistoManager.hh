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
// $Id: HistoManager.hh,v 1.13 2011-01-06 15:56:39 vnivanch Exp $
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

  void bookHisto();

  void BeginOfRun();
  void EndOfRun(G4int runID);

  void BeginOfEvent();
  void EndOfEvent();

  void ScoreNewTrack(const G4Track* aTrack);
  void AddEnergy(G4double edep, G4int idx, G4int copyNo);

  void AddDeltaElectron(const G4DynamicParticle*);
  void AddPhoton(const G4DynamicParticle*);

  G4double GetTrackLength() const {return trackLength;};
  void ResetTrackLength() {trackLength = 0.0, trackAbs = true;};
  void SetTrackOutAbsorber() {trackAbs = false;};
  G4bool GetTrackInAbsorber() const {return trackAbs;};
  void AddTrackLength(G4double x)   {trackLength += x;};
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
  void AddStep() {n_step += 1.0;};

  // Acceptance parameters
  void SetEdepAndRMS(G4int, G4ThreeVector);
  void SetBeamEnergy(G4double val) {beamEnergy = val;};

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
  G4double n_step;
  G4double n_step_target;
  G4bool trackAbs;        // Track is in absorber
  G4int n_evt;
  G4int n_elec;
  G4int n_posit;
  G4int n_gam;
  G4int n_gamph;
  G4int n_gam_tar;
  G4int n_lowe;
  G4int nBinsE, nBinsEA, nBinsED;
  G4bool nTuple;

  G4double Eabs1, Eabs2, Eabs3, Eabs4;
  G4double     E[25];
  G4DataVector Evertex;
  G4DataVector Nvertex;
  G4DataVector brem;
  G4DataVector phot;
  G4DataVector comp;
  G4DataVector conv;

  G4double  edeptrue[3];
  G4double  rmstrue[3];
  G4double  limittrue[3];
  G4double  edep[6];
  G4double  erms[6];
  G4double  edeptr[6];
  G4double  ermstr[6];
  G4int     stat[6];
  G4int     nmax;

  Histo*    histo;
  


};

#endif
