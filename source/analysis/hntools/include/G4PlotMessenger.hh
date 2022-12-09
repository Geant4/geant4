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

// The messenger class for batch plotting.
// It implements commands in /analysis/plot directory.
//
// Author: Ivana Hrivnacova, 21/10/2015  (ivana@ipno.in2p3.fr)

#ifndef G4PlotMessenger_h
#define G4PlotMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>
#include <string_view>

class G4PlotParameters;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIparameter;

class G4PlotMessenger : public G4UImessenger
{
  public:
    explicit G4PlotMessenger(G4PlotParameters* plotParameters);
    G4PlotMessenger() = delete;
    ~G4PlotMessenger() override;

    // Methods
    void SetNewValue(G4UIcommand* command, G4String value) final;

  private:
    // Helper functions
    template <typename CMD>
    std::unique_ptr<CMD> CreateCommand(G4String name, G4String guidance);
    void AddIntParameter(
      G4UIcommand& command, G4String name, G4String guidance, G4String range = "");

    // Functions to create commands
    void SetStyleCmd();
    void SetLayoutCmd();
    void SetDimensionsCmd();

    // constants
    static constexpr std::string_view fkClass { "G4PlotMessenger" };

    // Data members
    G4PlotParameters*  fPlotParameters { nullptr }; ///< Associated class
    std::unique_ptr<G4UIdirectory>  fDirectory;

    std::unique_ptr<G4UIcommand>  fSetLayoutCmd;
    std::unique_ptr<G4UIcommand>  fSetDimensionsCmd;
    std::unique_ptr<G4UIcmdWithAString>  fSetStyleCmd;
};

//_____________________________________________________________________________
template <typename CMD>
std::unique_ptr<CMD> G4PlotMessenger::CreateCommand(
  G4String name, G4String guidance)
{
  G4String fullName = "/analysis/plot/" + name;

  auto command = std::make_unique<CMD>(fullName, this);
  command->SetGuidance(guidance.c_str());
  command->AvailableForStates(G4State_PreInit, G4State_Idle);

  return command;
}

#endif

