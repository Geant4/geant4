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
#ifndef G4DNAEventScheduler_hh
#define G4DNAEventScheduler_hh 1

#include <G4VScheduler.hh>
#include <vector>
#include <map>
#include <memory>
#include "globals.hh"
#include "G4ITModelHandler.hh"
#include "G4ITStepStatus.hh"
#include "G4ITTrackHolder.hh"
#include "G4VStateDependent.hh"
#include "G4ITReaction.hh"
#include "G4DNAEventSet.hh"
#include "G4DNAUpdateSystemModel.hh"
#include "G4H2O2.hh"
#include <CLHEP/Units/SystemOfUnits.h>
class G4VITStepModel;
class G4DNAGillespieDirectMethod;
class G4UserMeshAction;
class IEventScheduler
{
 public:
  IEventScheduler()          = default;
  virtual ~IEventScheduler() = default;
};

class G4DNAEventScheduler : public IEventScheduler
{
 public:
  using MolType    = const G4MolecularConfiguration*;
  using MapList    = std::map<MolType, size_t>;
  using MapCounter = std::map<MolType, G4int>;
  G4DNAEventScheduler(const G4DNABoundingBox& boundingBox, G4int pixel);
  ~G4DNAEventScheduler() override;
  G4DNAEventScheduler(const G4DNAEventScheduler&) = delete;
  G4DNAEventScheduler& operator=(const G4DNAEventScheduler& right) = delete;
  void Initialize();
  void InitializeInMesh();
  void Voxelizing();
  void ReVoxelizing(G4int);
  void SetEndTime(const G4double&);
  G4double GetStartTime() const;
  G4double GetEndTime() const;
  G4double GetTimeStep() const;
  void SetStartTime(G4double time);

  inline void SetVerbose(G4int verbose) { fVerbose = verbose; }
  inline G4int GetVerbose() const;
  void Stepping();
  void SetChangeMesh(G4bool change) { fSetChangeMesh = change; }

  void Reset();
  void ResetInMesh();
  void RunInMesh();
  void Run();

  void AddTimeToRecord(const G4double& time);

  void RecordTime();
  void ClearAndReChargeCounter();
  void PrintRecordTime();
  void Stop();
  void SetMaxNbSteps(G4int);
  std::map<G4double /*time*/, MapCounter> GetCounterMap() const;
  G4DNAMesh* GetMesh() const;
  G4int GetPixels() const;
  void SetUserMeshAction(std::unique_ptr<G4UserMeshAction>);
  static G4bool CheckingReactionRadius(G4double resolution);

 private:
  G4int fVerbose = 0;
  G4bool fInitialized = false;
  G4double fStartTime = 1 * CLHEP::picosecond;
  G4double fEndTime = 10000 * CLHEP::second;
  G4int fStepNumber = 0;
  G4int fMaxStep = INT_MAX;
  G4bool fRunning = true;
  G4double fTimeStep = DBL_MAX;
  G4double fGlobalTime = 1 * CLHEP::picosecond;
  G4double fJumpingNumber = 0;
  G4double fReactionNumber = 0;
  G4int fPixel;
  G4bool fIsChangeMesh = false;
  G4bool fSetChangeMesh = true;
  G4int fStepNumberInMesh = 0;
  G4double fInitialPixels;
  G4double fTransferTime = 0.;
  const G4double C = 20;
  const G4double D = G4H2O2::Definition()->GetDiffusionCoefficient();  // this is the biggest D

  std::unique_ptr<G4DNAMesh> fpMesh;
  std::unique_ptr<G4DNAGillespieDirectMethod> fpGillespieReaction;
  std::unique_ptr<G4DNAEventSet> fpEventSet;
  std::unique_ptr<G4DNAUpdateSystemModel> fpUpdateSystem;
  std::unique_ptr<G4UserMeshAction> fpUserMeshAction;
  // an aternative Counter
  std::map<G4double /*time*/, MapCounter> fCounterMap;
  std::set<G4double> fTimeToRecord;
  std::set<G4double>::iterator fLastRecoredTime;
};
#endif
