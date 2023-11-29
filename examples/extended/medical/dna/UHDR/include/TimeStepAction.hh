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
#ifndef TimeStepAction_h
#define TimeStepAction_h

#include "G4DNAMesh.hh"
#include "G4UserTimeStepAction.hh"
#include <set>

class G4DNAEventScheduler;

class G4VPrimitiveScorer;

class PulseAction;

class G4VChemistryWorld;

class TimeStepAction : public G4UserTimeStepAction {
public:
  explicit TimeStepAction(const G4VChemistryWorld*, PulseAction* pPulse = nullptr);

  ~TimeStepAction() override = default;

  TimeStepAction(const TimeStepAction &other) = delete;

  TimeStepAction &operator=(const TimeStepAction &other) = delete;

  void StartProcessing() override { ; }

  void UserPreTimeStepAction() override;

  void UserPostTimeStepAction() override;

  void UserReactionAction(const G4Track & /*trackA*/,
                          const G4Track & /*trackB*/,
                          const std::vector<G4Track *> * /*products*/) override;

  void EndProcessing() override;

  G4DNAEventScheduler *GetEventScheduler() const;

  void SetInitialPixel();

private:
  std::unique_ptr<G4DNAEventScheduler> fpEventScheduler;

  void CompartmentBased();

  PulseAction* fpPulse = nullptr;
  const G4VChemistryWorld* fpChemWorld = nullptr;
  G4int fPixel = 0;

};

#endif
