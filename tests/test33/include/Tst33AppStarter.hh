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
//
// $Id: Tst33AppStarter.hh,v 1.12 2008-04-21 09:00:03 ahoward Exp $
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
#include "G4PlaceOfAction.hh"

class Tst33PhysicsList;
class Tst33DetectorConstruction;
class Tst33VGeometry;
class Tst33ParallelGeometry;
class G4VSampler;
//class G4GeometrySampler;
class G4VIStore;
class G4VWeightWindowStore;
class G4CellStoreScorer;
class G4CellScorer;
class G4CellScorerStore;
class G4UserRunAction;
class Tst33VEventAction;
class G4ProcessPlacer;
class Tst33WeightChangeProcess;
class G4VWeightWindowAlgorithm;

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
  void CreateWeightWindowStore(G4PlaceOfAction poa,
			       G4bool zeroWindow);
  void CreateWeightRoulette(G4int mode);
  void ClearSampling();
  void Run(G4int nevents);
  void AddWeightChanger();

  inline void ForcingCoupled(G4bool vl=true)
   { forceCoupled = vl; }

private:
  Tst33AppStarter(const Tst33AppStarter &);
  Tst33AppStarter &operator=(const Tst33AppStarter &);

  G4bool IsGeo_n_App();
  G4bool CheckCreateScorer();
  void DeleteScorers();
  G4bool CheckCreateIStore();
  G4bool CheckCreateWeightWindowStore();
  G4bool CheckCreateWeightRoulette();
  G4bool CheckCreateApp();
  void NewFailed(const G4String &function, const G4String &cl);

  Tst33AppStarterMessenger fMessenger;
  G4RunManager fRunManager;
  Tst33VGeometry *fMassGeometry;
  Tst33VGeometry *fSampleGeometry;
  //  Tst33VGeometry *fParallelGeometry;
  Tst33ParallelGeometry *fParallelGeometry;
  Tst33VApplication *fApp;
  Tst33DetectorConstruction *fDetectorConstruction;
  G4VSampler *fSampler;
  //  G4GeometrySampler *fSampler;
  G4CellStoreScorer *fScorer;
  G4CellScorerStore *fScorerStore;
  const G4CellScorer *fCell_19_Scorer;
  G4VIStore *fIStore;
  G4VWeightWindowStore *fWWStore;
  Tst33VEventAction *fEventAction;
  G4UserRunAction *fRunAction;
  G4bool fConfigured;
  G4bool fWeightroulette;
  G4int fTime;
  G4ProcessPlacer *fChangeWeightPlacer;
  Tst33WeightChangeProcess *fWeightChangeProcess;
  G4VWeightWindowAlgorithm *fWWAlg;

  Tst33PhysicsList * physlist;

  G4bool parallel_geometry;

  G4bool forceCoupled;

};

#endif
