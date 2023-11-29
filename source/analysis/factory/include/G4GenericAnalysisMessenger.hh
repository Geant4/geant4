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

// The messenger class for G4GenericAnalysisManager.
// It implements command:
// - /analysis/setDefaultFileType
//
// Author: Ivana Hrivnacova, 18/10/2022  (ivana@ipno.in2p3.fr)

#ifndef G4GenericAnalysisMessenger_h
#define G4GenericAnalysisMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class G4GenericAnalysisManager;
class G4UIcmdWithAString;

class G4GenericAnalysisMessenger : public G4UImessenger
{
  public:
    explicit G4GenericAnalysisMessenger(G4GenericAnalysisManager* manager);
    G4GenericAnalysisMessenger() = delete;
    ~G4GenericAnalysisMessenger() override;

    // Methods
    void SetNewValue(G4UIcommand* command, G4String value) final;

  private:
    template <typename CMD>
    std::unique_ptr<CMD> CreateCommand(
      G4String name, G4String guidance, G4String paremeterName,
      G4bool ommitable = false);

    // Data members
    G4GenericAnalysisManager* fManager { nullptr }; ///< Associated class
    std::unique_ptr<G4UIcmdWithAString>  fSetDefaultFileTypeCmd;
};

//_____________________________________________________________________________
template <typename CMD>
std::unique_ptr<CMD> G4GenericAnalysisMessenger::CreateCommand(
  G4String name, G4String guidance, G4String paremeterName, G4bool ommitable)
{
  G4String fullName = "/analysis/" + name;

  auto command = std::make_unique<CMD>(fullName, this);
  command->SetGuidance(guidance.c_str());
  command->SetParameterName(paremeterName.c_str(), ommitable);
  command->AvailableForStates(G4State_PreInit, G4State_Idle);

  return command;
}

#endif
