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

// The messenger class for ntuple management.
// It implements commands in /analysis/ntuple directory.
// It is asscoiciated with G4VAnalysisManager and this delegates
// call to both ntuple booking and ntuple managers.
//
// Author: Ivana Hrivnacova, 05/05/2015  (ivana@ipno.in2p3.fr)

#ifndef G4NtupleMessenger_h
#define G4NtupleMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>
#include <string_view>

class G4VAnalysisManager;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class G4NtupleMessenger : public G4UImessenger
{
  public:
    explicit G4NtupleMessenger(G4VAnalysisManager* manager);
    G4NtupleMessenger() = delete;
    ~G4NtupleMessenger() override;

    // Methods
    void SetNewValue(G4UIcommand* command, G4String value) final;

  private:
    // Methods
    template <typename CMD>
    std::unique_ptr<CMD> CreateCommand(G4String name, G4String guidance);
    void AddIdParameter(G4UIcommand& command);

    void SetActivationCmd();
    void SetActivationToAllCmd();
    void SetFileNameCmd();
    void SetFileNameToAllCmd();
    void ListCmd();

    // Static data members
    static constexpr std::string_view fkClass { "G4NtupleMessenger" };

    // Data members
    G4VAnalysisManager*  fManager { nullptr }; ///< Associated class

    std::unique_ptr<G4UIdirectory>      fNtupleDir;
    std::unique_ptr<G4UIcommand>        fSetActivationCmd;
    std::unique_ptr<G4UIcmdWithABool>   fSetActivationAllCmd;
    std::unique_ptr<G4UIcommand>        fSetFileNameCmd;
    std::unique_ptr<G4UIcmdWithAString> fSetFileNameAllCmd;
    std::unique_ptr<G4UIcommand>        fListCmd;
};

//_____________________________________________________________________________
template <typename CMD>
std::unique_ptr<CMD> G4NtupleMessenger::CreateCommand(
  G4String name, G4String guidance)
{
  G4String fullName = "/analysis/ntuple/" + name;

  auto command = std::make_unique<CMD>(fullName, this);
  command->SetGuidance(guidance.c_str());
  command->AvailableForStates(G4State_PreInit, G4State_Idle);

  return command;
}

#endif
