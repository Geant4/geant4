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
//
// $Id: Tst33AppStarter.hh,v 1.4 2002-11-20 13:09:15 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33AppStarter
//
// Class description:
//
// Different options for the application are supported:
// mass, or parallel geometry, using visualization or timing
// and configureable sampling.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33AppStarter_hh
#define Tst33AppStarter_hh Tst33AppStarter_hh

#include "globals.hh"
#include "Tst33AppStarterMessenger.hh"
#include "G4RunManager.hh"

class Tst33DetectorConstruction;
class Tst33VGeometry;
class G4VSampler;
class G4VIStore;
class G4CellStoreScorer;
class G4CellScorer;
class G4CellScorerStore;
class G4UserRunAction;
class Tst33VEventAction;

class Tst33AppStarter {
public:
  Tst33AppStarter();
  ~Tst33AppStarter();
  void CreateVisApplication();
  void CreateTimedApplication(G4int timed);
  void CreateMassGeometry();
  void CreateParallelGeometry();
  void ConfigureSampling();
  void PostRun();
  void CreateScorer();
  void CreateIStore();
  void CreateWeightRoulette();
  void ClearSampling();
  void Run(G4int nevents);


private:
  Tst33AppStarter(const Tst33AppStarter &);
  Tst33AppStarter &operator=(const Tst33AppStarter &);

  G4bool IsGeo_n_App();
  G4bool CheckCreateScorer();
  void DeleteScorers();
  G4bool CheckCreateIStore();
  G4bool CheckCreateWeightRoulette();
  G4bool CheckCreateApp();
  void NewFailed(const G4String &function, const G4String &cl);

  Tst33AppStarterMessenger fMessenger;
  G4RunManager fRunManager;
  Tst33VGeometry *fMassGeometry;
  Tst33VGeometry *fSampleGeometry;
  Tst33VGeometry *fParallelGeometry;
  Tst33VApplication *fApp;
  Tst33DetectorConstruction *fDetectorConstruction;
  G4VSampler *fSampler;
  G4CellStoreScorer *fScorer;
  G4CellScorerStore *fScorerStore;
  const G4CellScorer *fCell_19_Scorer;
  G4VIStore *fIStore;
  Tst33VEventAction *fEventAction;
  G4UserRunAction *fRunAction;
  G4bool fConfigured;
  G4bool fWeightroulette;
  G4int fTime;
};

#endif
