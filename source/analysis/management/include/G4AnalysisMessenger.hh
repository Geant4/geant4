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

// The messenger class for G4VAnalysisManager.
// It implements commands:
// - /analysis/closeFile
// - /analysis/compression
// - /analysis/list
// - /analysis/openFile
// - /analysis/reset
// - /analysis/setActivation
// - /analysis/setFileName
// - /analysis/setHistoDirName
// - /analysis/setNtupleDirName
// - /analysis/verbose
// - /analysis/write
//
// Author: Ivana Hrivnacova, 24/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisMessenger_h
#define G4AnalysisMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class G4VAnalysisManager;
class G4NtupleMessenger;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class G4AnalysisMessenger : public G4UImessenger
{
  public:
    explicit G4AnalysisMessenger(G4VAnalysisManager* manager);
    G4AnalysisMessenger() = delete;
    ~G4AnalysisMessenger() override;

    // Methods
    void SetNewValue(G4UIcommand* command, G4String value) final;

  private:
    template <typename CMD>
    std::unique_ptr<CMD> CreateCommand(
      G4String name, G4String guidance, G4String paremeterName,
      G4bool ommitable = false);
    std::unique_ptr<G4UIcmdWithoutParameter> CreateCommandWithoutParameter(
      G4String name, G4String guidance);

    // Data members
    G4VAnalysisManager* fManager { nullptr }; ///< Associated class
    std::unique_ptr<G4NtupleMessenger>  fNtupleMessenger;

    std::unique_ptr<G4UIdirectory>         fAnalysisDir;
    std::unique_ptr<G4UIcmdWithAString>    fOpenFileCmd;
    std::unique_ptr<G4UIcmdWithoutParameter>  fWriteCmd;
    std::unique_ptr<G4UIcmdWithoutParameter>  fResetCmd;
         // These command need to revert order of execution master-worker
         // To be investigated
    std::unique_ptr<G4UIcmdWithABool>      fCloseFileCmd;
    std::unique_ptr<G4UIcmdWithABool>      fListCmd;
    std::unique_ptr<G4UIcmdWithABool>      fSetActivationCmd;
    std::unique_ptr<G4UIcmdWithAnInteger>  fVerboseCmd;
    std::unique_ptr<G4UIcmdWithAnInteger>  fCompressionCmd;
    std::unique_ptr<G4UIcmdWithAString>    fSetFileNameCmd;
    std::unique_ptr<G4UIcmdWithAString>    fSetHistoDirNameCmd;
    std::unique_ptr<G4UIcmdWithAString>    fSetNtupleDirNameCmd;
};

//_____________________________________________________________________________
template <typename CMD>
std::unique_ptr<CMD> G4AnalysisMessenger::CreateCommand(
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
