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
// Author: Christian Velten (2025)

#ifndef G4MOLECULECOUNTERMANAGERMESSENGER_HH
#define G4MOLECULECOUNTERMANAGERMESSENGER_HH 1

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

#include <memory>

class G4MoleculeCounterManager;

class G4MoleculeCounterManagerMessenger : public G4UImessenger
{
  public:
    explicit G4MoleculeCounterManagerMessenger(G4MoleculeCounterManager*);
    ~G4MoleculeCounterManagerMessenger() override = default;
    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    void InitializeCommands();

  private:
    G4MoleculeCounterManager* fpManager{nullptr};
    std::unique_ptr<G4UIdirectory> fpManagerDir;
    std::unique_ptr<G4UIcmdWithABool> fpActiveCmd;
    std::unique_ptr<G4UIcmdWithABool> fpResetBeforeEventCmd;
    std::unique_ptr<G4UIcmdWithABool> fpResetBeforeRunCmd;
    std::unique_ptr<G4UIcmdWithABool> fpAccumulateIntoMasterCmd;
    std::unique_ptr<G4UIcmdWithAnInteger> fpVerboseCmd;
};

#endif  // G4MOLECULECOUNTERMANAGERMESSENGER_HH
