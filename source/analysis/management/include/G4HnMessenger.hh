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

// The messenger class for histogram information management.
// It implements commands in /analysis/h1 directory.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4HnMessenger_h
#define G4HnMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class G4HnManager;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class G4HnMessenger : public G4UImessenger
{
  public:
    explicit G4HnMessenger(G4HnManager& manager);
    G4HnMessenger() = delete;
    ~G4HnMessenger() override;

    // Methods
    void SetNewValue(G4UIcommand* command, G4String value) final;

  private:
    // Helper functions
    G4String GetObjectType() const;
    template <typename CMD>
    std::unique_ptr<CMD> CreateCommand(G4String name, G4String guidance);
    void AddIdParameter(G4UIcommand& command);
    void AddOptionParameter(G4UIcommand& command, G4String optionName);

    void SetHnAsciiCmd();
    void SetHnActivationCmd();
    void SetHnActivationToAllCmd();
    void SetHnPlottingCmd();
    void SetHnPlottingToAllCmd();
    void SetHnFileNameCmd();
    void SetHnFileNameToAllCmd();

    // constants
    static constexpr std::string_view fkClass { "G4HnMessenger" };

    // Data members
    G4HnManager& fManager; ///< Associated class
    G4String fHnType;
    std::unique_ptr<G4UIcommand>        fSetAsciiCmd;
    std::unique_ptr<G4UIcommand>        fSetActivationCmd;
    std::unique_ptr<G4UIcmdWithABool>   fSetActivationAllCmd;
    std::unique_ptr<G4UIcommand>        fSetPlottingCmd;
    std::unique_ptr<G4UIcmdWithABool>   fSetPlottingAllCmd;
    std::unique_ptr<G4UIcommand>        fSetFileNameCmd;
    std::unique_ptr<G4UIcmdWithAString> fSetFileNameAllCmd;
};

//_____________________________________________________________________________
template <typename CMD>
std::unique_ptr<CMD> G4HnMessenger::CreateCommand(
  G4String name, G4String guidance)
{
  G4String fullName = "/analysis/" + fHnType + "/" + name;
  G4String fullGuidance = guidance + GetObjectType();

  auto command = std::make_unique<CMD>(fullName, this);
  command->SetGuidance(fullGuidance);
  command->AvailableForStates(G4State_PreInit, G4State_Idle);

  return command;
}

#endif
