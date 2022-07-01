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
#ifndef MOLECULAR_TIMESTEP_ACTION_HH
#define MOLECULAR_TIMESTEP_ACTION_HH

#include "G4UserTimeStepAction.hh"

class DNAGeometry;

class EventAction;

class G4ITTrackHolder;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TimeStepAction : public G4UserTimeStepAction
{
 public:
  explicit TimeStepAction(EventAction*);

  ~TimeStepAction() override;

  void StartProcessing() override;

  void UserReactionAction(const G4Track&, const G4Track&,
                          const std::vector<G4Track*>*) override;

  void UserPreTimeStepAction() override;

  void UserPostTimeStepAction() override;

  void RadicalKillDistance();

 protected:
  EventAction* fEventAction;
  DNAGeometry* fDNAGeometry = nullptr;
  G4double fRadicalKillDistance;
  G4ITTrackHolder* fpChemistryTrackHolder;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_TIMESTEP_ACTION_HH
