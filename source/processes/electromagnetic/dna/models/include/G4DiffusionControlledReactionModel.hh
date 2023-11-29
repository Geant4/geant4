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
// Author: Hoang TRAN : 21/2/2019

#ifndef G4DiffusionControlledReactionModel_hh
#define G4DiffusionControlledReactionModel_hh 1

#include "G4VDNAReactionModel.hh"
#include <vector>
class G4DNAMolecularReactionData;
class G4DiffusionControlledReactionModel : public G4VDNAReactionModel
{
 public:
  G4DiffusionControlledReactionModel();
  ~G4DiffusionControlledReactionModel() override;

  G4DiffusionControlledReactionModel(
    const G4DiffusionControlledReactionModel&) = delete;
  G4DiffusionControlledReactionModel& operator =(
    const G4DiffusionControlledReactionModel&) = delete;

  void Initialise(const G4MolecularConfiguration*, const G4Track&) override;
  void InitialiseToPrint(const G4MolecularConfiguration*) override;
  G4double GetReactionRadius(const G4MolecularConfiguration*,
                             const G4MolecularConfiguration*) override;
  G4double GetReactionRadius(const G4int&) override;

  G4bool FindReaction(const G4Track&, const G4Track&,
                      G4double /*reactionRadius*/,
                      G4double& /*separationDistance*/,
                      G4bool /*alongStepInteraction*/) override
  {
    return true;
  }
  G4double GetTimeToEncounter(const G4Track& trackA, const G4Track& trackB);

 private:
  const std::vector<const G4DNAMolecularReactionData*>* fpReactionData =
    nullptr;
};
#endif