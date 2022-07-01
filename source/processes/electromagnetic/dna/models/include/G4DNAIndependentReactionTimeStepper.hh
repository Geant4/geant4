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
// 20/2/2019
// Author: HoangTRAN

#ifndef G4DNAIndependentReactionTimeStepper_hh
#define G4DNAIndependentReactionTimeStepper_hh 1

#include "G4VITTimeStepComputer.hh"
#include "G4KDTreeResult.hh"
#include "G4IRTUtils.hh"
#include <memory>
#include <set>
#include "G4ITTrackHolder.hh"
#include "G4ITReaction.hh"
#include "G4ReferenceCast.hh"

class G4VDNAReactionModel;
class G4DNAMolecularReactionTable;
class G4MolecularConfiguration;
class G4Molecule;
class G4ITReactionSet;
class G4ITReactionChange;
class G4VITReactionProcess;
class G4ITTrackHolder;

class G4DNAIndependentReactionTimeStepper : public G4VITTimeStepComputer
{
 public:
  G4DNAIndependentReactionTimeStepper();
  ~G4DNAIndependentReactionTimeStepper() override = default;
  G4DNAIndependentReactionTimeStepper(
    const G4DNAIndependentReactionTimeStepper&) = delete;
  G4DNAIndependentReactionTimeStepper& operator =(
    const G4DNAIndependentReactionTimeStepper&) = delete;

  void Prepare() override;
  G4double CalculateStep(const G4Track&, const G4double&) override;
  G4double CalculateMinTimeStep(G4double, G4double) override;

  void SetReactionModel(G4VDNAReactionModel*);
  G4VDNAReactionModel* GetReactionModel();

  std::unique_ptr<G4ITReactionChange> FindReaction(
    G4ITReactionSet* pReactionSet, const G4double& currentStepTime = 0,
    const G4double& previousStepTime       = 0,
    const G4bool& reachedUserStepTimeLimit = false);
  void SetReactionProcess(G4VITReactionProcess* pReactionProcess);
  void SetVerbose(G4int);

 private:
  void InitializeForNewTrack();
  class Utils;
  void CheckAndRecordResults(const Utils& utils);

  G4double GetTimeToEncounter(const G4Track& trackA, const G4Track& trackB);

  G4bool fHasAlreadyReachedNullTime = false;
  const G4DNAMolecularReactionTable*& fMolecularReactionTable =
    reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable);
  G4VDNAReactionModel* fReactionModel     = nullptr;
  G4ITTrackHolder* fpTrackContainer       = G4ITTrackHolder::Instance();
  G4ITReactionSet* fReactionSet           = G4ITReactionSet::Instance();
  G4int fVerbose                          = 0;
  G4double fRCutOff                       = G4IRTUtils::GetRCutOff();
  G4VITReactionProcess* fpReactionProcess = nullptr;
  std::map<G4int, G4ThreeVector> fSampledPositions;
  std::set<G4int> fCheckedTracks;

  class Utils
  {
   public:
    Utils(const G4Track& tA, const G4Track& tB);
    ~Utils() = default;
    const G4Track& fTrackA;
    const G4Track& fTrackB;
    const G4Molecule* fpMoleculeA;
    const G4Molecule* fpMoleculeB;
  };
};
#endif
